"""
CIDD Generation Module - Core molecular generation functionality.

This module implements the main Generation class that orchestrates
LLM-guided molecular design using structure-based drug design principles.
"""

import os
import sys
import json
import time
import random
import subprocess
import requests
from typing import Optional, Tuple, List, Any
from time import sleep

from rdkit import Chem
from rdkit.Chem import AllChem

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from generation.vina_eval import vina_dock_crossdocked as vina_dock_custom
from generation.preprocess import load_sdf_and_add_polar_hs
from generation.get_fragment import frag_mol_brics
from utils.unimap_ret import ret_fragments
from agent_utils import Agent, Part


from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex


def generate_interaction_txt(protein_pdb_path):


    # Initialize PLIP analysis
    my_plip = PDBComplex()

    # Load the protein PDB file
    #print("555", protein_pdb_path)
    my_plip.load_pdb(protein_pdb_path)


    # Perform the analysis
    my_plip.analyze()

    # Generate and print the report
    #print("666", my_plip.interaction_sets)
    res = {}
    for key in sorted(my_plip.interaction_sets):
        binding_site = my_plip.interaction_sets[key]
        report = BindingSiteReport(binding_site)
        report_txt = report.generate_txt()
        
        res[key] = ";".join(report_txt)
    
    # convert to plain text

   
    return res

try:
    from volcenginesdkarkruntime import Ark
    HAS_VOLCENGINE = True
except ImportError:
    HAS_VOLCENGINE = False


def iupac_2_smiles(iupac_name: str) -> Optional[str]:
    """
    Convert IUPAC name to SMILES using STOUT API.

    Args:
        iupac_name: IUPAC chemical name

    Returns:
        SMILES string if successful, None otherwise
    """
    try:
        encoded_iupac_name = requests.utils.quote(iupac_name)
        api_url = f"https://stout.api.decimer.ai/latest/stout/IUPAC2SMILES?input_text={encoded_iupac_name}&converter=opsin&visualize=2D"

        response = requests.get(api_url, timeout=10)
        response.raise_for_status()

        data = response.json()
        return data.get("SMILES")

    except (requests.RequestException, KeyError, json.JSONDecodeError):
        return None


def smiles_2_iupac(smiles: str) -> Optional[str]:
    """
    Convert SMILES to IUPAC name using STOUT API.

    Args:
        smiles: SMILES string

    Returns:
        IUPAC name if successful, None otherwise
    """
    try:
        api_url = "https://stout.api.decimer.ai/latest/stout/SMILE2IUPAC?retranslate=false&format=text"
        headers = {
            "accept": "application/json",
            "Content-Type": "text/plain",
        }

        response = requests.post(api_url, headers=headers, data=smiles, timeout=10)
        response.raise_for_status()

        iupac_name = response.text.strip()
        return iupac_name if iupac_name else None

    except requests.RequestException:
        return None


def iupac_2_smiles_api(iupac_name):
    encoded_iupac_name = requests.utils.quote(iupac_name)
    api_url = f"https://stout.api.decimer.ai/latest/stout/IUPAC2SMILES?input_text={encoded_iupac_name}&converter=opsin&visualize=2D"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data.get("SMILES").startswith("Failed"):
                return None
            return data.get("SMILES")
        return None
    except requests.RequestException as e:
        print(f"API failed for {iupac_name}: {e}")
        return 0

def iupac_2_smiles_cli(iupac_name, opsin_path="/home/gaobowen/opsin-cli-2.8.0-jar-with-dependencies.jar"):
    try:
        process = subprocess.Popen(
            ["java", "-jar", opsin_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=iupac_name)
        if process.returncode == 0:
            if stdout.strip()!="":
                return stdout.strip()
            else:
                return None
        print(f"CLI error for {iupac_name}: {stderr}")
        return None
    except Exception as e:
        print(f"CLI failed for {iupac_name}: {e}")
        return None

def iupac_2_smiles_all(iupac_name, fallback_to_cli=True):
    smiles = iupac_2_smiles_api(iupac_name)
    if smiles==0 and fallback_to_cli:
        print(f"Falling back to CLI for {iupac_name}")
        smiles = iupac_2_smiles_cli(iupac_name)
    return smiles


try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False


def send_request(context: List[dict], message: str, temperature: float = 0.7, n: int = 1) -> Optional[str]:
    """
    Send a request to LLM API with proper error handling.
    
    Supports two modes:
    1. Azure OpenAI (Default): Set AZURE_OPENAI_API_KEY, AZURE_OPENAI_URL
    2. DeepSeek (Volcengine Ark): Set VOLCENGINE_ARK_API_KEY (and optionally VOLCENGINE_MODEL)

    Args:
        context: List of message dictionaries for conversation context
        message: The message to send
        temperature: Temperature parameter for response generation
        n: Number of completions to generate

    Returns:
        Response content from the API, or None if failed
    """
    context.append({"role": "user", "content": message})

    try:
        # Mode 1: Azure OpenAI (Default, using requests)
        if os.environ.get("AZURE_OPENAI_API_KEY"):
            api_key = os.environ.get("AZURE_OPENAI_API_KEY")
            url = os.environ.get("AZURE_OPENAI_URL")
            
            if not url:
                raise ValueError("AZURE_OPENAI_URL must be set")
            
            headers = {
                'api-key': api_key,
                'Content-Type': 'application/json'
            }
            
            payload = {"messages": context, "n": n, "temperature": temperature}
            
            try:
                response = requests.post(url, json=payload, headers=headers, timeout=180)
            except:
                return None
            
            if response.status_code == 200:
                response_data = response.json()
                content = response_data['choices'][0]['message']['content']
                context.append({"role": "assistant", "content": content})
                print(content)
                
                if "usage" in response_data:
                    usage = response_data["usage"]
                    prompt_tokens = usage.get("prompt_tokens", 0)
                    completion_tokens = usage.get("completion_tokens", 0)
                    total_tokens = usage.get("total_tokens", 0)
                    print(f"Tokens used — prompt: {prompt_tokens}, completion: {completion_tokens}, total: {total_tokens}")
                else:
                    print("Token usage info not available.")
                
                return content
            else:
                print(response.text)
                print(message)
                return None

        # Mode 2: DeepSeek via Volcengine Ark
        elif os.environ.get("VOLCENGINE_ARK_API_KEY"):
            if not HAS_VOLCENGINE:
                raise ImportError("volcenginesdkarkruntime package not available. Install with: pip install volcengine-python-sdk[ark]")
            
            api_key = os.environ.get("VOLCENGINE_ARK_API_KEY")
            model = os.environ.get("VOLCENGINE_MODEL", "deepseek-v3-241226")
            
            client = Ark(api_key=api_key)
            
            response = client.chat.completions.create(
                model=model,
                messages=context,
                temperature=temperature,
                timeout=180
            )
            
            content = response.choices[0].message.content
            if content:
                context.append({"role": "assistant", "content": content})
                print(content)
                return content
            return None

        else:
            raise ValueError("No API key found. Set AZURE_OPENAI_API_KEY or VOLCENGINE_ARK_API_KEY")

    except Exception as e:
        print(f"API request failed: {e}")
        return None


class InteractionAgent(Agent):
    """
    Agent responsible for analyzing protein-ligand interactions.

    This agent performs molecular docking and interaction analysis
    to evaluate binding affinity and interaction patterns.
    """

    def __init__(self, name: str, data_root: str = None):
        """
        Initialize the InteractionAgent.

        Args:
            name: Name of the agent
            data_root: Root directory for protein/ligand data files
        """
        super().__init__(name)
        self.data_root = data_root or os.environ.get('CIDD_DATA_ROOT', './data')

    def combine_pdb_files(self, protein_file: str, ligand_file: str, output_file: str):
        """
        Combine protein and ligand PDB files into a single file.

        Args:
            protein_file: Path to protein PDB file
            ligand_file: Path to ligand PDB file
            output_file: Path to output combined PDB file
        """
        try:
            with open(output_file, 'w') as out:
                for pdb_file in [protein_file, ligand_file]:
                    with open(pdb_file, 'r') as f:
                        out.write(f.read())
        except IOError as e:
            print(f"Error combining PDB files: {e}")

    def analyze_inter(self, ligand_filename: str, smi: str, idx: int, ligand: bool = False):
        """
        Analyze protein-ligand interactions.

        Args:
            ligand_filename: Path to reference ligand file
            smi: SMILES string of molecule to analyze
            idx: Index for naming temporary files
            ligand: Whether to include ligand in analysis

        Returns:
            Tuple of (interaction_analysis, docking_score, interaction_text, mol_object)
        """
        protein_fn = os.path.join(
            os.path.dirname(ligand_filename),
            os.path.basename(ligand_filename)[:10] + '.pdb'
        )

        protein_root = os.path.join(self.data_root, os.path.dirname(ligand_filename))
        protein_path = os.path.join(self.data_root, protein_fn)
        ligand_path = os.path.join(self.data_root, ligand_filename)

        thread_name = str(idx)

        try:
            score = vina_dock_custom(
                smi, protein_root, protein_path, ligand_path,
                thread_name, f"{protein_root}/{thread_name}.sdf",
                exhaustiveness=16
            )
            print(f"Docking score: {score}")
        except Exception as e:
            print(f"Docking failed: {e}")
            return None, None, None, None

        input_file = f"{protein_root}/{thread_name}.sdf"
        output_file = f"{protein_root}/{thread_name}.pdb"

        try:
            command = ["obabel", input_file, "-O", output_file]
            subprocess.run(command, check=True, capture_output=True)
            print(f"Format conversion successful for {thread_name}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Format conversion failed: {e}")
            return None, None, None, None

        try:
            mol = Chem.MolFromMolFile(input_file)
            if mol is None:
                print("Failed to load molecule from SDF")
                return None, None, None, None
        except Exception as e:
            print(f"Error loading molecule: {e}")
            return None, None, None, None

        frags, labels, frag_text = frag_mol_brics(mol)

        if ligand:
            self.combine_pdb_files(protein_path, f"{protein_root}/ligand.pdb", f"{protein_root}/agent_combine_{thread_name}.pdb")
        else:
            self.combine_pdb_files(
                protein_path,
                f"{protein_root}/{thread_name}.pdb",
                f"{protein_root}/agent_combine_{thread_name}.pdb"
            )

        combined_pdb = f"{protein_root}/agent_combine_{thread_name}.pdb"
        interaction_res = generate_interaction_txt(combined_pdb)
        
        interaction_text = ""
        for key in interaction_res:
            if key.startswith("UNL"):
                interaction_text += f"{key}\n{interaction_res[key]}\n"

        content = f"""The interaction analysis results are as follows:
{interaction_text}

The atoms, coordinates and labels of the fragments are as follows:
{frag_text}

The SMILES string of the molecule is: {smi}

[Important] Using the information from the interaction analysis and the fragment labels, based on the coordinates matching, you need to identify the key fragment on the molecule that conduct the interaction with the protein pocket.

[Important] The format should be: which fragment can interact with which residue in the protein pocket. For example: The fragment 'ethanol' can interact with the residue 'ASP123' in the protein pocket.

[Important] It should be at fragment level, instead of atom level. You need to use the information from the fragment labels to identify which atom belongs to which fragment.
"""
        
        context = []
        frag_inter_match = send_request(context, content)

        return score, interaction_text, frag_inter_match, mol

    def execute(self, ligand_filename, smi, idx, ligand=False):
        """Execute interaction analysis and get LLM evaluation."""
        score, interaction_text, frag_inter_match, mol = self.analyze_inter(
            ligand_filename, smi, idx, ligand
        )

        if frag_inter_match is None:
            return (
                "no response",
                "unknown because no conformation was generated and there was no successful binding",
                "Failure to generate conformations",
                None
            )

        content = f"""The interaction analysis results are as follows:
{interaction_text}

The information on which fragment can interact with which residue in the protein pocket is as follows:
{frag_inter_match}

The SMILES string of the molecule is: {smi}
The docking score of the molecule is: {score}

Based on the interaction, give analysis on the molecule-protein interaction. Please use critical thinking to analyze, pointing out both the good and the bad points.

[Important] Your response should also analyze what are key fragments in the molecule that are important for the interaction with the protein pocket.
"""

        response = send_request(self.context, content)
        return response, score, interaction_text, mol


cur_dir = os.path.dirname(os.path.abspath(__file__))

class DesignAgent(Agent):
    def __init__(self, name):
        super().__init__(name)
        with open(f"{cur_dir}/prior_knowledge.txt", "r") as f:
            self.prior_knowledge = f.read()
        with open(f"{cur_dir}/scaffold_hopping.txt", "r") as f:
            self.scaffold_hopping = f.read()
        with open(f"{cur_dir}/side_chain_methods.txt", "r") as f:
            self.side_chain = f.read()
        with open(f"{cur_dir}/experience_summary.txt", "r") as f:
            self.experience_summary = f.read()
        self.trajectory = []
        self.reflections = []
        self.analysereports = []
        self.designs = []
        self.scores = []
        self.first_report = ""
        self.first_score = 0
        self.context = []
        self.protein_analysis = None
        self.ret_frags = None

        self.icl_experience = None

    def pocket_analysis_agent(self, pocket_structure):


        prompt = "Given a pocket structure from the pdb file, please think what kind of molecules can bind to this pocket. You need to think the properties needed for the molecules. You can also provide some analysis of the pocket."
        prompt += "[Important] Your analysis will be used by another agent to design the molecules. \n"
        prompt += "[Requirements] Your analysis should be based on the exact pocket structure. Do not use hypothetical structures or general information. \n"
        prompt += f"The pocket structure is: \n {pocket_structure}. \n"
        prompt += "Please think step by step and provide your analysis."

        context = []
        response = send_request(context, prompt)

        return response


    def icl_summarize(self, protein_name, idx):

        prev_experience = ""

        path = f"/drug/gaobowen/cidd_generation_results/{protein_name}"

        count = 0

        prev_dirs = os.listdir(path)
        # shuffle with seed 2025

        random.seed(2025)
        random.shuffle(prev_dirs)


        for prev_dir in prev_dirs:
            if count>0:
                break
            files = os.listdir(os.path.join(path, prev_dir))
            if f"molecule_{idx}.txt" in files:
                prev_experience += "#" * 20 + "\n"
                prev_experience += f"Experience {count}: \n"
                with open(os.path.join(path, prev_dir, f"molecule_{idx}.txt"), "r") as f:
                    prev_experience += f.read()
                    f.close()
                prev_experience + "\n"*5
                count += 1

        prompt = "Here are some previous designs of this molecule. We provide you with the previous designs, interaction reports and docking scores. Based on the previous designs, generate insights and provide suggestions for the next design. \n"
        prompt += "This is the in context learning task. You need to summarize the experience and provide suggestions for the next design. \n"
        prompt += "The previous designs and reflections are as follows: \n" + prev_experience + "\n"

        context = []

        print(prompt)
        response = send_request(context, prompt)


        return response


    def experience_summarize(self, smis, reports):
        content = "The previous designs and reflections are as follows: \n"
        for i in range(len(smis)):
            content += "The smiles of the molecule is: " + smis[i] + "\n" \
                    + "The reflection on the design is: \n" + reports[i] + "\n"
        content += "Based on the previous designs and reflections, summarize the experience and provide suggestions for the next design. \n"

        response = send_request([], content)
        if response is None:
            return "Failure to generate experience summary"
        return response


    def fragment_anaylsis(self, smi):

        mol = Chem.MolFromSmiles(smi)

        frags, labels, frag_text = frag_mol_brics(mol)

        dic = {}

        for frag in frags:
            frag_mol = Chem.MolFromSmiles(frag)
            if frag_mol is not None:

                ring_info = frag_mol.GetRingInfo()
                ring_count = ring_info.NumRings()

                if ring_count > 0:
                    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
                    avg_ring_size = sum(ring_sizes) / ring_count
                else:
                    avg_ring_size = 0
                ring_count = str(int(ring_count))
                avg_ring_size = str(int(avg_ring_size))
                if f"ring_num_{ring_count}" in self.drugbank_frags:
                    if f"avg_ring_size_{avg_ring_size}" in self.drugbank_frags[f"ring_num_{ring_count}"]:
                        ret_frags = self.drugbank_frags[f"ring_num_{ring_count}"][f"avg_ring_size_{avg_ring_size}"]
                        dic[frag] = ret_frags

        res_dic = ret_fragments(dic)

        print(res_dic)


        return res_dic


    def design(self, protein_name, idx, protein_structure, ori_smi, ori_report, prev_smis, prev_reports, change_scaffold=False):


        # if self.protein_analysis == None:
        #     self.protein_analysis = self.pocket_analysis_agent(protein_structure)
        ret_frags = None

        # if self.ret_frags == None:
        #     ret_frags = self.fragment_anaylsis(ori_smi)
        #     self.ret_frags = ret_frags

        #print(ret_frags)

        #ret_frags = None

        context = []

        # Build the design prompt
        content_parts = [
            "[Instruction]",
            f"This is the original molecule: {ori_smi}",
            f"This is the interaction analysis: {ori_report}",
            "Based on the pocket analysis and interaction analysis, design modifications to the original molecule.",
            ""
        ]

        if change_scaffold:
            content_parts.append(
                "Your job is scaffold hopping: change the core scaffold structure while retaining key pharmacophores."
            )

            if ret_frags is not None:
                content_parts.append(
                    "There are also some druglike fragments that can be used as references for your design."
                )
                content_parts.append(
                    "Here are the original fragments in the molecule and a list of similar druglike fragments:"
                )
                for frag in ret_frags:
                    content_parts.append(f"Original fragment: {frag}")
                    content_parts.append(f"Similar druglike fragments: {ret_frags[frag]}")
                content_parts.append(
                    "You don't have to directly use the druglike fragments, but you can use them as references."
                )
        else:
            content_parts.extend([
                "Your job is side chain modification: change the type or position of substituent groups.",
                "Do not change its whole molecular structure or properties.",
                "For example: substituting a carboxyl group with a tetrazole group, replacing 3-F with 2-F, or replacing 3-F with 3-Cl.",
                f"The side chain modification should be based on the following methods: {self.side_chain}"
            ])

        content_parts.extend([
            "",
            "[Design Objectives]",
            "1. The modified molecule should be more stable and easier to synthesize. New fragments should be common and stable.",
            "2. Retain key properties of the original molecule (shape, size, functionality).",
            "3. Maintain key interaction types with the same residues as the original molecule.",
            "4. Replace uncommon or unstable fragments with more common and stable counterparts (e.g., replace cyclohexadiene with benzene).",
            "5. The modified molecule should be more drug-like than the original.",
            "6. Must not contain aniline, hydrazine/azo, or Michael acceptor groups. Toxic structural alerts should be optimized out.",
            "",
            "[Design Requirements]",
            "1. Do not modify the molecule too much.",
            "2. Only give modification steps - do not generate SMILES yourself. Another agent will generate the new molecule based on your design.",
            ""
        ])

        if len(prev_reports) > 0:
            experience = self.experience_summarize(prev_smis, prev_reports)
            content_parts.extend([
                "[Experience Summary]",
                "We have done some modifications before. Learn from previous designs and reflections.",
                f"Summary of previous designs and reflections: {experience}",
                ""
            ])

        if self.icl_experience is not None:
            content_parts.extend([
                "[In Context Learning Experience]",
                self.icl_experience,
                ""
            ])

        content_parts.append("Now, give us the modification suggestions and requirements:")
        content = "\n".join(content_parts)
        response = send_request(context, content)


        return response

    def reflection(self, protein_structure, smi, score, design, gen_report, ori_smi, ori_report, ori_score):
        """Generate reflection on the molecular modification."""
        content = f"""[Instruction]
===============================================================================
This is the original molecule structure: {ori_smi}
This is the interaction analysis of this molecule and the protein pocket:
{ori_report}
===============================================================================
This is the modification design generated by you: {design}
This is the molecule structure designed by another agent based on your modification: {smi}
This is the interaction analysis generated by interaction expert agent:
{gen_report}

Based on the interaction analysis, reflect on the modification and provide feedback on the design.
Previously, we asked you to modify a molecule to get a new molecule that binds better to the protein pocket.
We analyzed the binding of the new molecule to the protein pocket and compared it with the binding before modification.
Based on these analyses, you need to reflect on the previous modification and provide feedback on the design.
Besides, we hope you can find out which fragments inside the structure of these two molecules are useful for binding to the protein pocket.
===============================================================================
Please use critical thinking to analyze, pointing out both the good and the bad points. Your reflection is:"""


        response = send_request([], content)


        return response

    def choose(self, protein_structure, smiles, scores, reports):
        content = "[Instruction]\n"
        content += "Given a protein pocket and a number of different molecules, I need you to find a molecule that is best suited to bind this protein pocket.  "\
                + "We will provide you with interaction reports for each molecule and protein pocket and ask you to make choices based on this information. " \
                + "You should consider both the binding analysis and whether the molecule has potential to be a real drug. \n" \
                + "You should consider that the designed molecules must not contain aniline, hydrazine/azo, or Michael acceptor groups. Toxic structural alerts should be identified and optimized out. \n" \
                + "You should output the reason and the smiles string of the chosen molecule. The smiles string should be enclosed in a pair of $$$. For example: ...(The reason) The chosen smiles is: $$$CCO$$$\n\n"

        for i in range(len(reports)):
            content += "This is the molecule structure: \n" + smiles[i] + "\n" \
                    + "This is the interaction analysis of this molecule and the protein pocket: \n" + reports[i] + "\n"
            content += "===============================================================================\n\n"


        content += "Please think step by step: "


        chosen_smi = ""
        failed_count = -1
        while chosen_smi not in smiles:
            failed_count+=1
            if failed_count > 5:
                break
            context = []
            response = send_request(context, content)
            try:
                chosen_smi = response.split("$$$")[1]
            except:
                continue
        if chosen_smi not in smiles:

            # random select index
            index = random.randint(0, len(smiles)-1)
            chosen_smi = smiles[index]
            chosen_report = reports[index]
            chosen_score = scores[index]
            return chosen_smi, chosen_report, chosen_score, index, response


        index = 0
        for i in range(len(smiles)):
            if smiles[i] == chosen_smi:
                chosen_report = reports[i]
                chosen_score = scores[i]
                index = i
                break
        return chosen_smi, chosen_report, chosen_score, index, response

# agent that generates the molecule based on modifications designed by the design agent

class MolGenAgent(Agent):
    def __init__(self, name):
        super().__init__(name)
        self.trajectory = []
        self.context = []
    def iupac_2_smiles(self, iupac_name):
        content = "Your job is to convert the IUPAC name to the SMILES string. \n" \
                + "The IUPAC name is: " + iupac_name + "\n" \
                + "The returned smiles should be enclosed in a pair of $$$. For example: $$$CCO$$$\n"

        response = send_request(self.context, content)

        try:
            smiles = response.split("$$$")[1]
            return smiles
        except:
            return None

    def generate_mol(self, design, ori_smi):
        ori_content = "This is the original molecule: \n" + ori_smi + "\n" \
                + "This is the modification designed by the design expert: \n" + design + "\n" \
                + "Change the molecule structure based on the modification designed by the design expert.\n" \
                + "Importantly, the modified molecule should be valid and make sense in the context of medicinal chemistry.\n"


        context = []
        response = send_request(context, ori_content, 1, 1)


        content = "Now please only return the smiles string of the generated molecule. The smiles string should be enclosed in a pair of $$$. For example: $$$cccc$$$\n" \

        response = send_request(context, content, 1, 1)


        mol = None
        failed_count = -1
        new_smi = "#"
        while mol is None:
            failed_count+=1
            print("failed count: ", failed_count)
            if failed_count > 5:
                return None, None

            try:
                temp_smi = response.split("$$$")[1].strip()
            except:
                content = "Now please only return the smiles string of the generated molecule. The smiles string should be enclosed in a pair of $$$. For example: $$$cccc$$$\n"
                response = send_request(context, content, 1, 1)
                continue


            if temp_smi is None:
                mol = None
            else:
                mol = Chem.MolFromSmiles(temp_smi)
            if mol is not None:
                new_smi = temp_smi

            else:
                content = f"Your previous generated molecule cannot be converted to a valid rdkit mol. For your new generation, please generate a valid molecule.\n" \
                        + "Check whether there are unkekulized atoms, unclosed rings, explicit valence, or other issues. \n" \
                        + "The generated new smiles string should be enclosed in a pair of $$$. For example:" \
                        + "The smiles string of the molecule is: $$$CCO$$$\n"
                response = send_request(context, content, 1, 1)


        mol = Chem.MolFromSmiles(new_smi)
        if mol is None:
            return None, None

        return new_smi, failed_count

def write_to_txt(data_list: List[Any], file_name: str):
    """
    Write a list of data to a text file.

    Args:
        data_list: List of items to write
        file_name: Output file name
    """
    with open(file_name, 'w') as file:
        for item in data_list:
            file.write(str(item) + '\n')


class Generation(Part):
    """
    Main CIDD Generation class that orchestrates the molecular generation pipeline.

    This class combines multiple agents (interaction, design, generation) to perform
    iterative molecular optimization guided by LLMs and structure-based evaluation.
    """

    def __init__(self, name: str, save_dir: str, data_root: str = None):
        """
        Initialize the Generation system.

        Args:
            name: Name identifier for this generation instance
            save_dir: Directory to save results
            data_root: Root directory for protein/ligand data files
        """
        super().__init__(name)
        self.save_dir = save_dir
        self.data_root = data_root or os.environ.get('CIDD_DATA_ROOT', './data')

        # Initialize agents
        self.inter_agent = InteractionAgent("Interaction Agent", self.data_root)
        self.design_agent = DesignAgent("Design Agent")
        self.molgen_agent = MolGenAgent("Molecule Generation Agent")

        # Add agents to the part
        self.add_agent(self.inter_agent)
        self.add_agent(self.design_agent)
        self.add_agent(self.molgen_agent)

    def new_execute_all(self, ligand_filename: str, ligand_smi: str, idx: int,
                       num_scaffold_explore: int = 10, num_side_chain_explore: int = 3):
        """
        Execute the complete molecular generation and optimization pipeline.

        Args:
            ligand_filename: Path to reference ligand file
            ligand_smi: Starting molecule SMILES
            idx: Index for naming output files
            num_scaffold_explore: Number of scaffold modifications to explore
            num_side_chain_explore: Number of side chain modifications to explore

        Returns:
            Tuple of generation results or None if failed
        """
        inter_agent = InteractionAgent("Interaction Agent", self.data_root)
        designs = []
        smis = []
        scores = []
        reports = []
        reflections = []

        protein_fn = os.path.join(
            os.path.dirname(ligand_filename),
            os.path.basename(ligand_filename)[:10] + '.pdb'
        )

        pocket = os.path.join(self.data_root, protein_fn)

        analyze_report,first_score, ori_inter_text, ori_mol = inter_agent.execute(ligand_filename, ligand_smi, idx, ligand=False)
        if analyze_report is None:
            return None
        try:
            first_score = float(first_score)
        except:
            return None


        first_report = analyze_report
        initial_score = first_score

        with open(pocket, "r") as f:
            protein_structure = f.read()

        failed_counts = []
        protein_name = os.path.dirname(ligand_filename)
        gen_mol_list = []
        
        for i in range(5):
            try:
                design = self.design_agent.design(
                    protein_name, idx, protein_structure, ligand_smi,
                    first_report, smis, reflections, change_scaffold=True
                )
                if design is None:
                    continue

                new_smi, failed_count = self.molgen_agent.generate_mol(design, ligand_smi)
                if new_smi is None:
                    continue

                analyze_report, score, inter_text, gen_mol = self.inter_agent.execute(
                    ligand_filename, new_smi, idx
                )

                if isinstance(score, str):
                    continue
                if analyze_report is None:
                    continue


                try:
                    reflection = self.design_agent.reflection(protein_structure, new_smi, score, design, analyze_report, ligand_smi, first_report, initial_score)
                except:
                    reflection = "failed to generate reflection"
                if reflection == None:
                    reflection = "failed to generate reflection"
                smis.append(new_smi)
                designs.append(design)
                scores.append(score)
                reports.append(analyze_report)
                reflections.append(reflection)
                failed_counts.append(failed_count)
                gen_mol_list.append(gen_mol)


                file_name = f"{self.save_dir}/molecule_{idx}.txt"

                with open(file_name, 'a') as file:
                    # original ligand
                    file.write("#" * 20 + "Scaffold Round " + str(i+1) + "#" * 20 + "\n")
                    file.write("Original ligand: " + ligand_smi + "\n")
                    # new ligand
                    file.write("New ligand: " + new_smi + "\n")
                    # design
                    # ######## Design: ########
                    file.write("############# Design: #############\n")
                    file.write(design + "\n")
                    file.write("############# Original Docking Score: #############\n")
                    file.write(str(initial_score) + "\n")
                    file.write("############# Original Interaction report: #############\n")
                    file.write(ori_inter_text + "\n")
                    # write docking score
                    file.write("############# New Docking Score: #############\n")
                    file.write(str(score) + "\n")
                    # interaction report
                    file.write("############# New Interaction Report: #############\n")
                    file.write(inter_text + "\n")
                    # failed count
                    file.write("############# Failed Count: #############\n")
                    file.write(str(failed_count) + "\n")

                f.close()
            except:
                continue


        if len(smis) >= 3:

            try:
                res = self.design_agent.choose(protein_structure, smis, scores, reports)
                if res is None:
                    with open(f"{self.save_dir}/molecule_{idx}.txt", 'a') as file:
                        file.write("\n" + "#" * 20 + "fail to choose" + "#" * 20 + "\n")
                    file.close()
                    return None
                chosen_smi, chosen_report, chosen_score, index, reason =  res
                chosen_failed_count = failed_counts[index]

                chosen_gen_mol = gen_mol_list[index]

                with open(f"{self.save_dir}/molecule_{idx}.txt", 'a') as file:
                    file.write("\n" + "#" * 20 + "Scaffold Chosen Molecule part" + "#" * 20 + "\n")
                    file.write("The chosen molecule is: " + chosen_smi + "\n")
                    file.write("The docking score of the chosen molecule is: " + str(chosen_score) + "\n")
                    file.write("The reason for choosing the molecule is: " + reason + "\n")

                file.close()

                # save ori_mol and chosen gen_mol same file

                # with open(f"{self.save_dir}/molecule_{idx}_final_mols.pkl", 'wb') as f:
                #     pickle.dump([ori_mol, chosen_gen_mol], f)


                return chosen_smi, chosen_report, chosen_score, ligand_smi, initial_score, chosen_failed_count, ori_mol, chosen_gen_mol
            except:
                # get molecule with best score
                index = scores.index(min(scores))
                chosen_smi = smis[index]
                chosen_report = reports[index]
                chosen_score = scores[index]
                chosen_failed_count = failed_counts[index]
                chosen_gen_mol = gen_mol_list[index]
                return chosen_smi, chosen_report, chosen_score, ligand_smi, initial_score, chosen_failed_count, ori_mol, chosen_gen_mol
        else:
            return None


        

        #
def process(ligand_filename, ligand_smi, i, save_dir):
    generation_module = Generation("Generation Module", save_dir)
    res = generation_module.new_execute_all(ligand_filename, ligand_smi, i)
    if res is None:
        print("failed")
        return


    return ligand_filename


if __name__ == "__main__":
    """
    Example usage of CIDD generation system.
    This is for demonstration purposes only.
    """
    import torch
    from tqdm import tqdm

    target = "7ksi"
    cur_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    save_dir = os.path.join("./results", target, cur_time)
    os.makedirs(save_dir, exist_ok=True)

    # Load example data (adjust path as needed)
    data_path = os.environ.get("CIDD_EXAMPLE_DATA", "./data/example_molecules.pt")

    try:
        mol_data = torch.load(data_path)

        total_times = []
        for idx in range(min(10, len(mol_data))):
            start = time.time()

            mol = mol_data[idx]["mol"]
            ligand_smi = Chem.MolToSmiles(mol)
            ligand_filename = mol_data[idx]["ligand_filename"]

            print(f"Processing molecule {idx}: {ligand_smi}")

            process(ligand_filename, ligand_smi, idx, save_dir)

            elapsed = time.time() - start
            total_times.append(elapsed)
            print(f"Completed {idx} in {elapsed:.2f} seconds.")

        avg_time = sum(total_times) / len(total_times) if total_times else 0
        print(f"\nAverage time per task: {avg_time:.2f} seconds.")

    except FileNotFoundError:
        print(f"Example data not found at {data_path}")
    except Exception as e:
        print(f"Error during processing: {e}")

