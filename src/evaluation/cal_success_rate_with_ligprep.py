# import os
# import pandas as pd
#
#
# def calculate_region_ratio(data: pd.DataFrame, region: dict) -> dict:
#     """
#     Calculate the success ratio of molecules based on defined property regions.
#     A molecule is considered successful if at least one row in its group meets all region criteria.
#     """
#     # Initialize counters
#     success_count = 0
#     total_molecules = 0
#
#     # Group the data by 'molecule'
#     grouped = data.groupby("molecule")
#
#     for molecule, group in grouped:
#         total_molecules += 1
#         molecule_success = False
#
#         for _, row in group.iterrows():
#             try:
#                 # Check if the row meets all region criteria
#                 success = all(
#                     region[prop][0] <= row[prop] <= region[prop][1]
#                     for prop in region if prop in row and pd.notna(row[prop])
#                 )
#
#                 if success:
#                     molecule_success = True
#                     break  # Exit loop early if any row is successful
#
#             except Exception as e:
#                 print(f"Error processing molecule '{molecule}': {e}")
#
#         if molecule_success:
#             success_count += 1
#
#     # Calculate success and failure ratios
#     success_ratio = success_count / total_molecules if total_molecules > 0 else 0
#     failure_ratio = 1 - success_ratio
#     # print(data)
#     print('total_molecules', 'success_ratio', 'success_count')
#     print(total_molecules, success_ratio, success_count)
#
#     return {
#         "success_ratio": success_ratio,
#         "failure_ratio": failure_ratio
#     }
#
#
# def main(main_dir, log_time, tmp_dir):
#     print(os.path.basename(tmp_dir))
#
#     os.makedirs(tmp_dir, exist_ok=True)
#     os.chdir(tmp_dir)
#
#     try:
#         ori_data = pd.read_csv('ori_mols.CSV', on_bad_lines='skip')
#         gen_data = pd.read_csv('gen_mols.CSV', on_bad_lines='skip')
#     except FileNotFoundError as e:
#         print(f"Error: {e}")
#         return
#
#     # # Debug loaded data
#     # print("Original Data:")
#     # print(ori_data.head())
#     # print("Generated Data:")
#     # print(gen_data.head())
#
#     # Rename columns if necessary
#     ori_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
#     gen_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
#
#     # Define the regions for each property
#     region = {
#         "#rtvFG": (0, 2),
#         "CIQPlogS": (-6.5, 0.5),
#         "QPlogS": (-6.5, 0.5),
#         "QPlogPo_w": (-2.0, 6.5),
#         "mol_MW": (130.0, 725.0),
#         "dipole": (1.0, 12.5),
#         "SASA": (300.0, 1000.0),
#         "FOSA": (0.0, 750.0),
#         "FISA": (7.0, 330.0),
#         "IP(eV)": (7.9, 10.5),
#         "EA(eV)": (-0.9, 1.7),
#         "#metab": (1, 8),
#         "PercentHumanOralAbsorption": (25, 100),
#         "PSA": (7.0, 200.0),
#         "RuleOfFive": (0, 3),
#     }
#
#     # Calculate success ratios directly on DataFrame
#     print('ori')
#     ori_region_ratios = calculate_region_ratio(ori_data, region)
#     print('gen')
#     gen_region_ratios = calculate_region_ratio(gen_data, region)
#
#
#     print("\n=== Region-Based Ratios ===")
#     print(f"Original data success ratio: {ori_region_ratios['success_ratio']:.2%}")
#     # print(f"Original data failure ratio: {ori_region_ratios['failure_ratio']:.2%}")
#     print(f"Generated data success ratio: {gen_region_ratios['success_ratio']:.2%}")
#     # print(f"Generated data failure ratio: {gen_region_ratios['failure_ratio']:.2%}")
#     print("="*30)
#     print()
#
#
# if __name__ == '__main__':
#     main(
#         main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
#         log_time='2025-01-01-12-03-52',
#         tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_bfn_cd_lp'
#     )
#
#
#     main(
#         main_dir='/mnt/nfs-ssd/data/gaobowen/td_pdbbind_0.6_pa_8',
#         log_time='2024-12-27-18-59-09',
#         tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_td_lp'
#     )
#
#     main(
#         main_dir='/mnt/nfs-ssd/data/gaobowen/p2m_pdbbind_0.6_pa_8',
#         log_time='2024-12-28-16-36-19',
#         tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_p2m_lp'
#     )
#
#
#
#
#
#
#

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


# def calculate_region_ratio(data: dict, region: dict) -> dict:
#     """Calculate the ratio of values falling inside and outside the region for each property."""
#     ratios = {}
#     for prop in data:
#         if prop not in region:
#             continue  # Skip properties without defined regions
#
#         try:
#             in_region = [x for x in data[prop] if region[prop][0] <= x <= region[prop][1]]
#             out_region = [x for x in data[prop] if x < region[prop][0] or x > region[prop][1]]
#
#             in_ratio = len(in_region) / len(data[prop]) if len(data[prop]) > 0 else 0
#             out_ratio = len(out_region) / len(data[prop]) if len(data[prop]) > 0 else 0
#
#             ratios[prop] = {
#                 "in_region": in_ratio,
#                 "out_region": out_ratio
#             }
#         except Exception as e:
#             print(f"Error processing property '{prop}': {e}")
#             ratios[prop] = {
#                 "in_region": 0,
#                 "out_region": 0
#             }
#
#     return ratios


def calculate_region_ratio(data: pd.DataFrame, region: dict) -> dict:
    """
    Calculate the success ratio of molecules based on defined property regions.
    A molecule is considered successful if at least one row in its group meets all region criteria.
    """
    # Initialize counters
    success_count = 0
    total_molecules = 0

    # Group the data by 'molecule'
    grouped = data.groupby("molecule")

    for molecule, group in grouped:
        total_molecules += 1
        molecule_success = False

        for _, row in group.iterrows():
            try:
                # Check if the row meets all region criteria
                success = all(
                    region[prop][0] <= row[prop] <= region[prop][1]
                    for prop in region if prop in row and pd.notna(row[prop])
                )

                if success:
                    molecule_success = True
                    break  # Exit loop early if any row is successful

            except Exception as e:
                print(f"Error processing molecule '{molecule}': {e}")

        if molecule_success:
            success_count += 1

    # Calculate success and failure ratios
    success_ratio = success_count / total_molecules if total_molecules > 0 else 0
    failure_ratio = 1 - success_ratio

    return {
        "success_ratio": success_ratio,
        "failure_ratio": failure_ratio
    }



# def select_best_fit_lines_within_region(df: pd.DataFrame, region: dict) -> pd.DataFrame:
#     """
#     For each molecule, select the line with the maximum number of properties within the defined regions.
#     """
#     best_fit_lines = []
#
#     # Group by molecule
#     grouped = df.groupby("molecule")
#     for molecule, group in grouped:
#         # Count properties within the region for each line
#         group["within_region_count"] = group.apply(
#             lambda row: sum(
#                 region[prop][0] <= row[prop] <= region[prop][1]
#                 for prop in region if prop in row and pd.notna(row[prop])
#             ),
#             axis=1,
#         )
#
#         # Select the line with the highest count of properties within the region
#         best_fit_line = group.loc[group["within_region_count"].idxmax()]
#         best_fit_lines.append(best_fit_line)
#
#     # Combine all best-fit lines into a new DataFrame
#     return pd.DataFrame(best_fit_lines)


# def calculate_success_rate(df: pd.DataFrame, region: dict) -> float:
#     """
#     Calculate the ratio of molecules where all specified properties fall within their regions.
#     A molecule is considered 'successful' if it meets all region criteria.
#     """
#     mask = pd.Series([True] * len(df))
#     for prop, (low, high) in region.items():
#         if prop in df.columns:
#             # Ensure the property is numeric
#             df[prop] = pd.to_numeric(df[prop], errors='coerce')
#             mask = mask & df[prop].between(low, high)
#         else:
#             mask = mask & False  # If property missing, set to False
#
#     success_count = mask.sum()
#     total = len(df)
#     ratio = success_count / total if total > 0 else 0
#     return ratio


def main(main_dir, log_time, tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_1'):
    os.makedirs(tmp_dir)
    os.chdir(tmp_dir)

    log_path = f'/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/meta_log_{log_time}.txt'
    print(f'Processing logs and saving new meta in {log_path}')

    ori_smis = []
    gen_smis = []

    with open(log_path, 'w') as out:
        for target in tqdm(os.listdir(main_dir), desc="Processing targets"):
            out.write(f'============= {target} =============\n')
            target_dir = os.path.join(main_dir, target, log_time)
            if os.path.isdir(target_dir):
                for f in os.listdir(target_dir):
                    if f.endswith('.pkl'):
                        pkl_path = os.path.join(target_dir, f)

                        with open(pkl_path, 'rb') as pf:
                            try:
                                chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, failed_count = pickle.load(
                                    pf)
                            except Exception as e:
                                print(f"Error loading pickle file '{pkl_path}': {e}")
                                continue

                            molecule_name = f.strip().split('.')[0]
                            out.write(f'{target}\t{molecule_name}\t{chosen_score}\t{ori_score}\n')

                            ori_smis.append(ori_smi)
                            gen_smis.append(chosen_smi)
            else:
                print(f'skipping {target} as it does not contain log_time {log_time}')
            out.write(f'\n')

    # ori_smi_path = os.path.join(tmp_dir, 'ori_mols.smi')
    # gen_smi_path = os.path.join(tmp_dir, 'gen_mols.smi')

    with open('ori_mols.smi', 'w') as ori_smi_f:
        for smi in ori_smis:
            ori_smi_f.write(f'{smi}\n')

    with open('gen_mols.smi', 'w') as gen_smi_f:
        for smi in gen_smis:
            gen_smi_f.write(f'{smi}\n')

    SCHRODINGER = '/opt/schrodinger2021-2'
    # os.chdir(tmp_dir)

    print(f'ligprep with {SCHRODINGER}')

    # lig_prep_dir = os.path.join(home, 'ligand_prepare_results')

    # os.chdir(lig_prep_dir)

    def lig_prep(SCHRODINGER, name, lig_conf, smi_path, conformation_out, n_conf=1, n_cpu=64, ):

        with open(lig_conf, 'w') as f:
            f.write(
                f'''INPUT_FILE_NAME {smi_path}
        OUT_MAE {conformation_out}
        FORCE FIELD 16
        EPIK yes
        DETERMINE CHIRALITIES yes
        RESPECT CHIRALITIES yes
        IGNORE CHIRALITIES no
        NUM STEREOISOMERS {n_conf}
        '''
            )
        # if not os.path.exists(conformation_out):
        os.system(
            f'{SCHRODINGER}/ligprep -inp {lig_conf} -HOST localhost:{n_cpu} -NJOBS {n_cpu} -JOBNAME ligprep_{name} -WAIT -LOCAL')

        # # convert maegz to sdf
        # os.system(
        #     f'{SCHRODINGER}/utilities/structconvert {os.path.splitext(smiles)[0] + ".maegz"} {os.path.splitext(smiles)[0] + ".sdf"}')

    # ori_mae_path = os.path.join(tmp_dir, 'ori_mols.maegz')
    # gen_mae_path = os.path.join(tmp_dir, 'gen_mols.maegz')

    lig_conf = os.path.join(tmp_dir, 'conf.in')

    #MARK

    lig_prep(SCHRODINGER, 'ori', lig_conf, 'ori_mols.smi', 'ori_mols.maegz')
    lig_prep(SCHRODINGER, 'gen', lig_conf, 'gen_mols.smi', 'gen_mols.maegz')

    print(f'Start calculating QikProp with {SCHRODINGER}')
    print('cwd:', os.getcwd())
    # Run QikProp for original molecules
    cmd_ori = f'{SCHRODINGER}/qikprop -i ori_mols.maegz'
    print(f'Running command: {cmd_ori}')
    os.system(cmd_ori)

    # Run QikProp for generated molecules
    cmd_gen = f'{SCHRODINGER}/qikprop -i gen_mols.maegz'
    print(f'Running command: {cmd_gen}')
    os.system(cmd_gen)

    print(f'Start parsing and plotting QikProp results...')

    # Define the regions for each property
    region = {
        "#rtvFG": (0, 2),
        "CIQPlogS": (-6.5, 0.5),
        "QPlogS": (-6.5, 0.5),
        "QPlogPo_w": (-2.0, 6.5),
        "mol_MW": (130.0, 725.0),
        "dipole": (1.0, 12.5),
        "SASA": (300.0, 1000.0),
        "FOSA": (0.0, 750.0),
        "FISA": (7.0, 330.0),
        # "QPPCaco": (25, 1000000),
        "IP(eV)": (7.9, 10.5),
        "EA(eV)": (-0.9, 1.7),
        "#metab": (1, 8),
        # "QPlogKhsa": (-1.5, 1.5),
        # "HumanOralAbsorption": (2, 3),
        "PercentHumanOralAbsorption": (25, 100),
        "PSA": (7.0, 200.0),
        "RuleOfFive": (0, 3),
    }

    # Properties to extract and compare (ensure all region keys are included)
    properties = [
        "#rtvFG", "QPlogPo_w", "QPlogS", "CIQPlogS",
        "QPlogBB", "#metab", "RuleOfFive", "#ringatoms",
        "mol_MW", "dipole", "SASA", "FOSA", "FISA",
        "QPPCaco", "IP(eV)", "EA(eV)",
        "QPlogKhsa", "HumanOralAbsorption",
        "PercentHumanOralAbsorption", "PSA"
    ]

    # Read QikProp output data
    try:
        ori_data = pd.read_csv('ori_mols.CSV', on_bad_lines='skip')
        gen_data = pd.read_csv('gen_mols.CSV', on_bad_lines='skip')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # Rename columns if necessary
    ori_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
    gen_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)

    # Check if data is not empty
    if ori_data.empty or gen_data.empty:
        print("Error: Failed to load QikProp data.")
        return

    # Extract properties
    ori_props = extract_properties(ori_data, properties)
    gen_props = extract_properties(gen_data, properties)

    # Plot comparisons (optional)
    # plot_comparisons(ori_props, gen_props, properties, fig_name_prefix="property_comparison")

    # Calculate the ratio of values inside and outside the region
    ori_region_ratios = calculate_region_ratio(ori_props, region)
    gen_region_ratios = calculate_region_ratio(gen_props, region)

    print("\n=== Region-Based Ratios ===")
    print("Original properties ratio in/out of regions:")
    for prop in ori_region_ratios:
        print(
            f"{prop}: In: {ori_region_ratios[prop]['in_region']:.2f}, Out: {ori_region_ratios[prop]['out_region']:.2f}")

    print("\nGenerated properties ratio in/out of regions:")
    for prop in gen_region_ratios:
        print(
            f"{prop}: In: {gen_region_ratios[prop]['in_region']:.2f}, Out: {gen_region_ratios[prop]['out_region']:.2f}")

    # # Apply this function before calculating success rates
    # ori_data_filtered = select_best_fit_lines_within_region(ori_data, region)
    # gen_data_filtered = select_best_fit_lines_within_region(gen_data, region)

    # # Calculate success rate (molecules passing all region requirements)
    # ori_success_rate = calculate_success_rate(ori_data_filtered, region)
    # gen_success_rate = calculate_success_rate(gen_data_filtered, region)

    # # ori_success_rate = calculate_success_rate(ori_data, region)
    # ori_success_rate = calculate_success_ratio(ori_data, region)
    # # gen_success_rate = calculate_success_rate(gen_data, region)
    # gen_success_rate = calculate_success_ratio(gen_data, region)
    #
    # print("\n=== Success Rate ===")
    # print(f"Original dataset success rate: {ori_success_rate:.2%}")
    # print(f"Generated dataset success rate: {gen_success_rate:.2%}")


if __name__ == '__main__':
    # # # Example of processing multiple log_times and tmp_dirs
    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-01-12-03-52',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_bfn_cd_lp'
    # )

    # # # # You can uncomment and add more `main` function calls as needed
    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/td_pdbbind_0.6_pa_8',
    #     log_time='2024-12-27-18-59-09',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_td_lp'
    # )

    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/p2m_pdbbind_0.6_pa_8',
    #     log_time='2024-12-28-16-36-19',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_p2m_lp'
    # )

    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-04-06-15-30',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_bfn_cd_lp_250104'
    # )


    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-05-16-19-34',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_td_lp_250106'
    # )

    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-07-13-24-54',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_decompdiff_250107'
    # )

    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-08-01-02-55',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_mc_4o_250108'
    # )

    main(
        main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
        log_time='2025-01-07-13-54-59',
        tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_p2m_mini_250107'
    )

    # main(
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-08-14-05-35',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_decompdiff_250108'
    # )