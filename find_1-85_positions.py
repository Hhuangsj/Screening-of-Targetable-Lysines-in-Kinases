import pandas as pd
def get_positions(excel_path):
    df = pd.read_excel(excel_path)
    data_dict = df.set_index('residue number').to_dict()['region']
    return data_dict
  
def find_acid(my_dict):
    acid_positions = {
        'D': {},
        'E': {},
        'H': {},
        'K': {},
        'M': {},
        'R': {},
        'S': {},
        'T': {},
        'W': {}
    }
    for i in range(1, 86):
        for acid_type in acid_positions:
            acid_positions[acid_type][i] = ""

    for name, seqs in my_dict.items():
        #print(seqs)
        for index, char in enumerate(seqs):
            if char == "D":
                acid_positions['D'][index + 1] += name+', '
            if char == "E":
                acid_positions['E'][index + 1] += name+', '
            if char == "H":
                acid_positions['H'][index + 1] += name+', '
            if char == "K":
                acid_positions['K'][index + 1] += name+', '
            if char == "M":
                acid_positions['M'][index + 1] += name+', '
            if char == "R":
                acid_positions['R'][index + 1] += name+', '
            if char == "S":
                acid_positions['S'][index + 1] += name+', '
            if char == "T":
                acid_positions['T'][index + 1] += name+', '
            if char == "W":
                acid_positions['W'][index + 1] += name+', '

                    
    return acid_positions

xlsx_file_path = "klifs_all.xlsx"  
date_file = pd.read_excel(xlsx_file_path)
klifs_dict = dict(zip(date_file.iloc[:, 0], date_file.iloc[:, 1]))
excel_path='/home/hsj/workstations/Analysis/kinase/positions/position_region.xlsx'
dict=get_positions(excel_path)
retules = find_acid(klifs_dict)
kineas_df = pd.DataFrame(retules)
file_path = "positions.xlsx"
kineas_df.to_excel(file_path, index=False)
