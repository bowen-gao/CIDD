import shutil
import os
import pickle
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import matplotlib.pyplot as plt


def extract_properties(df: pd.DataFrame, properties: list) -> dict:
    """Extract specified properties from the DataFrame."""
    extracted = {}
    for prop in properties:
        if prop in df.columns:
            extracted[prop] = df[prop].dropna().tolist()
        else:
            print(f"Warning: Property '{prop}' not found in the data.")
            extracted[prop] = []
    return extracted


def plot_comparisons(ori_props: dict, gen_props: dict, properties: list, fig_name_prefix: str = "comparison"):
    """Plot comparison figures for specified properties."""
    for prop in properties:
        if not ori_props[prop] or not gen_props[prop]:
            print(f"Skipping plot for '{prop}' due to missing data.")
            continue

        # Calculate mean values
        mean_ori = sum(ori_props[prop]) / len(ori_props[prop])
        mean_gen = sum(gen_props[prop]) / len(gen_props[prop])

        # Create the plot
        plt.figure(figsize=(8, 6))
        plt.hist(ori_props[prop], bins=40, alpha=0.5, label="Original", color="blue")
        plt.hist(gen_props[prop], bins=40, alpha=0.5, label="Generated", color="orange")
        plt.axvline(x=mean_ori, color="blue", linestyle="--", label=f"Mean (Original): {mean_ori:.2f}")
        plt.axvline(x=mean_gen, color="orange", linestyle="--", label=f"Mean (Generated): {mean_gen:.2f}")
        plt.title(f"Comparison of {prop}")
        plt.xlabel(f"{prop} Value")
        plt.ylabel("Frequency")
        plt.legend()
        plt.tight_layout()

        # Save the plot
        fig_name = f"{fig_name_prefix}_{prop}.png"
        plt.savefig(fig_name, dpi=300)
        print(f"Plot for '{prop}' saved as {fig_name}.")

        plt.close()


def main(main_dir, log_time, tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_1'):
    log_path = f'/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/meta_log_{log_time}.txt'
    print(f'Processing logs and save new meta in {log_path}')

    ori_smis = []
    gen_smis = []

    with open(log_path, 'w') as out:
        for target in tqdm(os.listdir(main_dir)):
            out.write(f'============= {target} =============\n')
            if log_time in os.listdir(os.path.join(main_dir, target)):
                dir = os.path.join(main_dir, target, log_time)
                for f in os.listdir(dir):
                    if f.endswith('.pkl'):
                        pkl_path = os.path.join(dir, f)

                        with open(pkl_path, 'rb') as pf:
                            # print(pickle.load(pf))
                            # chosen_iupac, chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, failed_count = pickle.load(pf)
                            chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, failed_count = pickle.load(pf)

                            molecule_name = f.strip().split('.')[0]
                            # print(molecule_name)
                            out.write(f'{target}\t{molecule_name}\t{chosen_score}\t{ori_score}\n')

                            ori_smis.append(ori_smi)
                            gen_smis.append(chosen_smi)
            else:
                print(f'skipping {target}')
            out.write(f'\n')

    print(f'Start generating .sdf files in {tmp_dir}')

    def write_sdf(sdf_path, smis_list):
        writer = Chem.SDWriter(sdf_path)
        for smi in tqdm(smis_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:  # If molecule creation is successful
                # Add hydrogens to the molecule
                mol = Chem.AddHs(mol)
                # Generate 3D coordinates
                try:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Use ETKDG method for 3D generation
                    # Optimize the conformation using MMFF force field
                    AllChem.MMFFOptimizeMolecule(mol)
                    writer.write(mol)
                except Exception as e:
                    print(f"Failed to generate 3D conformation for SMILES: {smi}. Error: {e}")
            else:
                print(f"Skipping invalid SMILES: {smi}")
        writer.close()
        print(f".sdf file has been saved to {sdf_path}")

    ori_sdf_path = os.path.join(tmp_dir, 'ori_mols.sdf')
    gen_sdf_path = os.path.join(tmp_dir, 'gen_mols.sdf')

    print(ori_sdf_path)
    write_sdf(ori_sdf_path, ori_smis)
    print(gen_sdf_path)
    write_sdf(gen_sdf_path, gen_smis)

    # ori_smi_path = os.path.join(tmp_dir, 'ori_mols.smi')
    # gen_smi_path = os.path.join(tmp_dir, 'gen_mols.smi')
    #
    # with open(ori_smi_path, 'w') as ori_smi_f:
    #     for smi in ori_smis:
    #         ori_smi_f.write(f'{smi}\n')
    #
    # with open(gen_smi_path, 'w') as gen_smi_f:
    #     for smi in gen_smis:
    #         gen_smi_f.write(f'{smi}\n')

    SCHRODINGER = '/opt/schrodinger2021-2'
    # print(f'ligprep with {SCHRODINGER}')

    # lig_prep_dir = os.path.join(home, 'ligand_prepare_results')
    # os.chdir(lig_prep_dir)

    # def lig_prep(SCHRODINGER, name, lig_conf, smi_path, conformation_out, n_conf=1, n_cpu=8, ):
    #
    #     with open(lig_conf, 'w') as f:
    #         f.write(
    #             f'''INPUT_FILE_NAME {smi_path}
    #     OUT_MAE {conformation_out}
    #     FORCE FIELD 16
    #     EPIK yes
    #     DETERMINE CHIRALITIES yes
    #     RESPECT CHIRALITIES yes
    #     IGNORE CHIRALITIES no
    #     NUM STEREOISOMERS {n_conf}
    #     '''
    #         )
    #     # if not os.path.exists(conformation_out):
    #     os.system(
    #         f'{SCHRODINGER}/ligprep -inp {lig_conf} -HOST localhost:{n_cpu} -NJOBS {n_cpu} -JOBNAME ligprep_{name} -WAIT -LOCAL')
    #
    #     # # convert maegz to sdf
    #     # os.system(
    #     #     f'{SCHRODINGER}/utilities/structconvert {os.path.splitext(smiles)[0] + ".maegz"} {os.path.splitext(smiles)[0] + ".sdf"}')
    #
    # ori_mae_path = os.path.join(tmp_dir, 'ori_mols.maegz')
    # gen_mae_path = os.path.join(tmp_dir, 'gen_mols.maegz')
    #
    # lig_conf = os.path.join(tmp_dir, 'conf.in')
    #
    # lig_prep(SCHRODINGER, 'ori', lig_conf, ori_smi_path, ori_mae_path)
    # lig_prep(SCHRODINGER, 'gen', lig_conf, gen_smi_path, gen_mae_path)

    print(f'Start calculating qikprop with {SCHRODINGER}')

    os.chdir(tmp_dir)
    print(os.getcwd())

    # os.system(f'{SCHRODINGER}/qikprop -i {ori_mae_path}')
    os.system(f'{SCHRODINGER}/qikprop -i ori_mols.sdf')
    # os.system(f'{SCHRODINGER}/qikprop -i {gen_mae_path}')
    os.system(f'{SCHRODINGER}/qikprop -i gen_mols.sdf')

    print(f'Start parsing and plotting qikprop...')

    # Properties to extract and compare
    properties = ["QPlogPo_w", "QPlogS", "CIQPlogS", "QPlogBB", "#metab", "RuleOfFive", "#ringatoms"]

    # Read data
    ori_data = pd.read_csv('ori_mols.CSV', on_bad_lines='skip')
    gen_data = pd.read_csv('gen_mols.CSV', on_bad_lines='skip')

    ori_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
    gen_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)

    # Check if data is not empty
    if ori_data.empty or gen_data.empty:
        print("Error: Failed to load data.")
        return

    # Extract properties
    ori_props = extract_properties(ori_data, properties)
    gen_props = extract_properties(gen_data, properties)

    # Plot comparisons
    plot_comparisons(ori_props, gen_props, properties, fig_name_prefix="property_comparison")


if __name__ == '__main__':
    main(
        main_dir='/mnt/nfs-ssd/data/gaobowen/bfn_pdbbind_0.6_pa_8',
        log_time='2024-12-26-22-27-57',
        tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir'
    )

