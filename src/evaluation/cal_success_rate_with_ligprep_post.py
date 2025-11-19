import os
import pandas as pd


def get_input_count(main_dir, log_time, tmp_dir):
    if tmp_dir == '/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_fda':
        return 2153
    else:
        input_count = 0
        for target in os.listdir(main_dir):
            target_dir = os.path.join(main_dir, target, log_time)
            if os.path.isdir(target_dir):
                for f in os.listdir(target_dir):
                    if f.endswith('.pkl'):
                        input_count += 1
            else:
                print(f'skipping {target} as it does not contain log_time {log_time}')
        return input_count


def calculate_region_ratio(data: pd.DataFrame, region: dict, main_dir, log_time, tmp_dir) -> dict:
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

    input_count = get_input_count(main_dir, log_time, tmp_dir)
    print(f'input count: {input_count}')
    with open(log_file_path, 'a') as log:
        log.write(f'        input count: {input_count}\n')

    # Calculate success and failure ratios
    success_ratio = success_count / input_count if input_count > 0 else 0
    failure_ratio = 1 - success_ratio
    # print(data)
    print('total_molecules', 'success_ratio', 'success_count')
    print(total_molecules, success_ratio, success_count)
    
    with open(log_file_path, 'a') as log:
        log.write('         total_molecules, success_ratio, success_count\n')
        log.write(f'        {total_molecules}, {success_ratio}, {success_count}\n')

    return {
        "success_ratio": success_ratio,
        "failure_ratio": failure_ratio
    }


def main(region, log_file_path, main_dir, log_time, tmp_dir):
    print(os.path.basename(tmp_dir))
    with open(log_file_path, 'a') as log:
        log.write(f"    {log_time}")
        log.write(f"    {os.path.basename(tmp_dir)}\n")

    os.makedirs(tmp_dir, exist_ok=True)
    os.chdir(tmp_dir)

    try:
        ori_data = pd.read_csv('ori_mols.CSV', on_bad_lines='skip')
        gen_data = pd.read_csv('gen_mols.CSV', on_bad_lines='skip')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # # # Debug loaded data
    # print("Original Data:")
    # print(ori_data.head())
    # print("Generated Data:")
    # print(gen_data.head())

    # Rename columns if necessary
    ori_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
    gen_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
    ori_data.rename(columns={"dip^2/V": "dip^2_V"}, inplace=True)
    gen_data.rename(columns={"dip^2/V": "dip^2_V"}, inplace=True)
    ori_data.rename(columns={"ACxDN^.5/SA": "ACxDN^.5_SA"}, inplace=True)
    gen_data.rename(columns={"ACxDN^.5/SA": "ACxDN^.5_SA"}, inplace=True)

    # # Define the regions for each property
    # region = {
    #     "#rtvFG": (0, 2),
    #     "CIQPlogS": (-6.5, 0.5),
    #     "QPlogS": (-6.5, 0.5),
    #     "QPlogPo_w": (-2.0, 6.5),
    #     "mol_MW": (130.0, 725.0),
    #     "dipole": (1.0, 12.5),
    #     "SASA": (300.0, 1000.0),
    #     "FOSA": (0.0, 750.0),
    #     "FISA": (7.0, 330.0),
    #     "IP(eV)": (7.9, 10.5),
    #     "EA(eV)": (-0.9, 1.7),
    #     "#metab": (1, 8),
    #     "PercentHumanOralAbsorption": (25, 100),
    #     "PSA": (7.0, 200.0),
    #     "RuleOfFive": (0, 3),
    # }

    # Calculate success ratios directly on DataFrame
    print('ori')
    with open(log_file_path, 'a') as log:
        log.write(f"    ori\n")
    ori_region_ratios = calculate_region_ratio(ori_data, region, main_dir, log_time, tmp_dir)
    print('gen')
    with open(log_file_path, 'a') as log:
        log.write(f"    gen\n")
    gen_region_ratios = calculate_region_ratio(gen_data, region, main_dir, log_time, tmp_dir)


    print("\n=== Region-Based Ratios ===")
    print(f"Original data success ratio: {ori_region_ratios['success_ratio']:.2%}")
    # print(f"Original data failure ratio: {ori_regio_n_ratios['failure_ratio']:.2%}")
    print(f"Generated data success ratio: {gen_region_ratios['success_ratio']:.2%}")
    # print(f"Generated data failure ratio: {gen_region_ratios['failure_ratio']:.2%}")
    print("="*30)
    print()
    with open(log_file_path, 'a') as log:
        log.write("         === Region-Based Ratios ===\n")
        log.write(f"        Original data success ratio: {ori_region_ratios['success_ratio']:.2%}\n")
        log.write(f"        Generated data success ratio: {gen_region_ratios['success_ratio']:.2%}\n")
        log.write(f"    " +"="*60+ '\n\n')


def fda_main(region, log_file_path, main_dir, log_time, tmp_dir):
    print(os.path.basename(tmp_dir))

    os.makedirs(tmp_dir, exist_ok=True)
    os.chdir(tmp_dir)

    try:
        ori_data = pd.read_csv('ori_mols.CSV', on_bad_lines='skip')
        # gen_data = pd.read_csv('gen_mols.CSV', on_bad_lines='skip')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # # # Debug loaded data
    # print("Original Data:")
    # print(ori_data.head())
    # print("Generated Data:")
    # print(gen_data.head())

    # Rename columns if necessary
    ori_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)
    ori_data.rename(columns={"dip^2/V": "dip^2_V"}, inplace=True)
    ori_data.rename(columns={"ACxDN^.5/SA": "ACxDN^.5_SA"}, inplace=True)
    # gen_data.rename(columns={"QPlogPo/w": "QPlogPo_w"}, inplace=True)

    # Calculate success ratios directly on DataFrame
    print('ori')
    ori_region_ratios = calculate_region_ratio(ori_data, region, main_dir, log_time, tmp_dir)
    # print('gen')
    # gen_region_ratios = calculate_region_ratio(gen_data, region)


    print("\n=== Region-Based Ratios ===")
    print(f"Original data success ratio: {ori_region_ratios['success_ratio']:.2%}")
    # print(f"Original data failure ratio: {ori_region_ratios['failure_ratio']:.2%}")
    # print(f"Generated data success ratio: {gen_region_ratios['success_ratio']:.2%}")
    # print(f"Generated data failure ratio: {gen_region_ratios['failure_ratio']:.2%}")
    print("="*30)
    print()
    
    with open(log_file_path, 'a') as log:
        log.write("         === Region-Based Ratios ===\n")
        log.write(f"        Original data success ratio: {ori_region_ratios['success_ratio']:.2%}\n")
        log.write(f"    " +"="*60+ '\n\n')

if __name__ == '__main__':
    # Define the regions for each property
    region2 = {
        # "#rtvFG": (0, 2),
        # # "CIQPlogS": (-6.5, 0.5),
        # "QPlogS": (-6.5, 0.5),
        # "QPlogPo_w": (-2.0, 6.5),
        # "mol_MW": (130.0, 725.0),
        # "dipole": (1.0, 12.5),
        # "SASA": (300.0, 1000.0),
        # "FOSA": (0.0, 750.0),
        # "FISA": (7.0, 330.0),
        # "IP(eV)": (7.9, 10.5),
        # "EA(eV)": (-0.9, 1.7),
        # "#metab": (1, 8),
        # "PercentHumanOralAbsorption": (25, 100),
        # "PSA": (7.0, 200.0),
        "RuleOfFive": (0, 0),
    }


    region1 = {  # origenal
        "#rtvFG": (0, 2),
        # "CIQPlogS": (-6.5, 0.5),
        "QPlogS": (-6.5, 0.5),
        "QPlogPo_w": (-2.0, 6.5),
        "mol_MW": (130.0, 725.0),
        "dipole": (1.0, 12.5),
        "SASA": (300.0, 1000.0),
        "FOSA": (0.0, 750.0),
        "FISA": (7.0, 330.0),
        "IP(eV)": (7.9, 10.5),
        "EA(eV)": (-0.9, 1.7),
        "#metab": (1, 8),
        "PercentHumanOralAbsorption": (25, 100),
        "PSA": (7.0, 200.0),
        # "RuleOfFive": (0, 3),
    }

    # Define the regions for each property
    region0 = {
        "#stars": (0, 5),
        "#amine": (0, 1),
        "#amidine": (0, 0),
        "#acid": (0, 1),
        "#amide": (0, 1),
        "#rotor": (0, 15),
        "#rtvFG": (0, 2),
        # "CNS": (-2, 2),
        "mol_MW": (130.0, 725.0),
        "dipole": (1.0, 12.5),
        "SASA": (300.0, 1000.0),
        "FOSA": (0.0, 750.0),
        "FISA": (7.0, 330.0),
        "PISA": (0.0, 450.0),
        "WPSA": (0.0, 175.0),
        "volume": (500.0, 2000.0),
        "donorHB": (0.0, 6.0),
        "accptHB": (2.0, 20.0),
        "dip^2/V": (0.0, 0.13),
        "ACxDN^.5/SA": (0.0, 0.05),
        "glob": (0.75, 0.95),
        "QPpolrz": (13.0, 70.0),
        "QPlogPC16": (4.0, 18.0),
        "QPlogPoct": (8.0, 35.0),
        "QPlogPw": (4.0, 45.0),
        "QPlogPo/w": (-2.0, 6.5),
        "QPlogS": (-6.5, 0.5),
        "CIQPlogS": (-6.5, 0.5),
        # "QPlogHERG": (-5, float("inf")),  # Concern below -5
        "QPPCaco": (25, float("inf")),  # Poor below 25
        "QPlogBB": (-3.0, 1.2),
        "QPPMDCK": (25, float("inf")),  # Poor below 25
        "QPlogKp": (-8.0, -1.0),
        "IP(eV)": (7.9, 10.5),
        "EA(eV)": (-0.9, 1.7),
        "#metab": (1, 8),
        "QPlogKhsa": (-1.5, 1.5),
        "PercentHumanOralAbsorption": (25, 100),  # High >80%
        "SAFluorine": (0.0, 100.0),
        "SAamideO": (0.0, 35.0),
        "PSA": (7.0, 200.0),
        "#NandO": (2, 15),
        # "RuleOfFive": (0, 4),
        "RuleOfThree": (0, 3),
        # "Jm": (0, float("inf")),  # No explicit upper limit given
    }

    # main(
    #     region,
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-01-12-03-52',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_bfn_cd_lp'
    # )
    #
    # main(
    #     region,
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/td_pdbbind_0.6_pa_8',
    #     log_time='2024-12-27-18-59-09',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_td_lp'
    # )
    #
    # main(
    #     region,
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/p2m_pdbbind_0.6_pa_8',
    #     log_time='2024-12-28-16-36-19',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_p2m_lp'
    # )

    # main(
    #     region,
    #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
    #     log_time='2025-01-04-06-15-30',
    #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_bfn_cd_lp_250104'
    # )

    log_file_path = '/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/codes/log_all_250108_1.txt'
    for i, region in enumerate([region0, region1, region2]):
        print(f"region = region{i}")
        with open(log_file_path, 'a') as log:
            log.write(f"region = region{i}\n\n")
        
        main(
            region,
            log_file_path,
            main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
            log_time='2025-01-05-16-19-34',
            tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_td_lp_250106'
        )
    
        # main(
        #     region,
        #     main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
        #     log_time='2025-01-07-13-24-54',
        #     tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_decompdiff_250107'
        # )
    
        main(
            region,
            log_file_path,
            main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
            log_time='2025-01-08-01-02-55',
            tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_mc_4o_250108',
            
        )
    
        main(
            region,
            log_file_path,
            main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
            log_time='2025-01-07-13-54-59',
            tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_p2m_mini_250107'
        )
    
        main(
            region,
            log_file_path,
            main_dir='/mnt/nfs-ssd/data/gaobowen/molcraft_crossdocked/',
            log_time='2025-01-08-14-05-35',
            tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_decompdiff_250108'
        )
    
        fda_main(
            region,
            log_file_path,
            main_dir='',
            log_time='',
            tmp_dir='/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/qikprop/tmp_dir_fda'
        )




