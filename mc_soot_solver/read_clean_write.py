# This file will read unformatted Cantera and FlameMaster output files and rewrite the data: 1 parameter to each row

import csv

distance = []
time = []
temp = []
conc = []

def read_clean_write(parameter, file_name):

    # Read in data
    f = open(file_name, 'r')
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        parameter = parameter + row

    # Check & delete '' values
    print(str(parameter))
    del parameter[-1]
    del parameter[-1]

    # Print again to compare
    print(parameter)

    return parameter

distance = read_clean_write(distance, file_name='ISF_flame_2b_distance_m.csv')
time = read_clean_write(time, file_name='ISF_flame_2b_time_ms.csv')
temp = read_clean_write(temp, file_name='ISF_flame_2b_temp_K.csv')
conc = read_clean_write(conc, file_name='ISF_flame_2b_A4_molefraction.csv')

with open('ISF_flame_2b_data.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerows([distance, time, temp, conc])
