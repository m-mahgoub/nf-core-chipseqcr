#!/usr/bin/env python

import sys
import yaml
import pandas as pd

with open(sys.argv[1], 'r') as file:
    yaml_dict = yaml.safe_load(file)


all_info_as_dictionary = {}
file_paths_list = []  # local and remote paths for all plots, not including files generated in the workflow

for plot, plot_meta in yaml_dict.items():
    plot_name = plot_meta['name']
    all_info_as_dictionary[plot_name] = {}
    for key, user_list in plot_meta.items():
        files_names_string = ''
        files_labels_string = ''
        if key == 'bed' or key == 'bigwig':
            for file_meta in user_list:
                file_label = file_meta[1]['label']
                file_source = list(file_meta[0].keys())[0]
                
                # bed file generated from the pipline
                if file_source == 'sampleID':  
                    file_name = file_meta[0][file_source] + '.' + key
                    if file_label == "":
                        file_label = file_meta[0][file_source]
                
                # bed file from local or remote path
                if file_source == 'file':
                    file_path = file_meta[0][file_source]
                    file_name = file_path.split('/')[-1]
                    if file_label == "":
                        file_label = file_name

                    if file_path not in file_paths_list:
                        file_paths_list.append(file_path)

                files_names_string += file_name + ' ' # add files names without quotes
                files_labels_string += f'"{file_label}" '  # add files names without quotes
            all_info_as_dictionary[plot_name][key + '_files_inLine'] = files_names_string
            all_info_as_dictionary[plot_name][key + '_labels_inLine_with_quotes'] = files_labels_string




df = pd.DataFrame.from_dict(all_info_as_dictionary, orient='index')
df.index.name = 'plot_name'

# save the inline strings to be injected in the commandline of deeptool process
df.to_csv(sys.argv[2],sep='\t',quotechar="'")

# save the paths for local and remote files used for deeptool process
df_paths = pd.DataFrame(file_paths_list)
df_paths.to_csv(sys.argv[3],sep='\t', index=False, header=False)
