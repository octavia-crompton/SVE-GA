# # coding=utf-8
# import json
# import os
# import pandas as pd

# def get_name_tuples(core):
#     names = []
#     for key in core.index:
#         sim = core.loc[key]
#         name = ", ".join([ "{0}={1}".format(name_var, sim[name_var]) for name_var in name_vars])
#         names.append((name, key))
#     return names

# def print_input_params(directory_name, include=None):
#     """

#     Parameters
#     ----------
#     directory_name
#     """
#     if include is None:
#         include = []
#     file_name = '{0}/all_params.json'.format(directory_name)
#     all_params = json.load(open(file_name))
#     print_all_params(all_params, include)


# def load_all_params(directory_name):
#     """
#     Parameters
#     ----------
#     directory_name
#     """
#     file_name = '{0}/all_params.json'.format(directory_name)
#     all_params = json.load(open(file_name))

#     return all_params


# def print_all_params(all_params, include=None):
#     """
#     Print factorial combinations of these variables
#     """
#     print('batch:')
#     for key in all_params['batch_dict']:
#         print('\t' + key + ' : ' + str(all_params['batch_dict'][key])[1:-1])
#     print('sim:')
#     for key in all_params['sim_dict']:
#         print('\t' + key + ' : ' + str(all_params['sim_dict'][key])[1:-1])
#     print('common:')
#     common_vars = ['H_i']
#     if include:
#         common_vars = common_vars + include
#     for key in common_vars:
#         if key in all_params['common_dict'].keys():
#             print('\t' + key + ' : ' + str(all_params['common_dict'][key]))


# def extract_params(sim):
#     """
#     Return sim parameters in a dictionary

#     Note:
#         DO not include flags
#     """
#     keys = ['q1_m2hr', 'Ks', 'dt_sw', 'tr', 'p', 'tmax_scale', 'dt_print',
#             'save_fluxes', 'save_sve', 'nrow', 'ncol', 'dx',
#             'veg_type', 'fV', 'grad_fV', 'seed', 'sigma_scale', 'sigma',
#             'stripe_count', 'downslope', 'spots',
#             'topo', 'So', 'imodel', 's_scale', 'theta_r', 'theta_s',
#             'theta_i', 'H_i', 'Ao', 'scheme', 'alpha', 'alphaB',
#             'itype1', 'itype3', 'itype2', 'itype4',  'epsh']
#     params = {}
#     for key in keys:
#         if key in sim.keys():
#             params[key] = sim[key]
#     return params


# def list_image_names(base_dir, core, remove_extension=True):
#     """

#     Parameters
#     ----------
#     base_dir
#     core
#     remove_extension

#     Returns
#     -------

#     image_names = list_image_names(base_dir, core)
#     """
#     image_list = os.listdir(os.path.join(base_dir, core.image_dir[0]))
#     if remove_extension:
#         image_list = [name.split('.')[0] for name in image_list]
#     return image_list


# def flatten_nested(nested_dict):
#     """
#     Flattens a nested dictionary

#     Parameters
#     ----------
#     nested_dict

#     Returns
#     -------

#     """
#     flattened_dict = {}
#     for _, item in nested_dict.items():
#         for key, nested_item in item.items():
#             if type(nested_item) == list:
#                 if len(nested_item) == 1:
#                     nested_item = nested_item[0]
#             flattened_dict[key] = nested_item

#     return pd.Series(flattened_dict)


# def summarize_param_files(project_dir):
#     """
#     Summarizes all the param files in the project directory Searches the
#     `model_output` subdirectory and summarizes all `all_params` files

#     Parameters
#     ----------
#     project_dir : path
#         Project directory

#     Returns
#     -------

#     """
#     output_dir = os.path.join(project_dir, 'model_output')

#     model_output_summary = pd.DataFrame()
#     for base_name in os.listdir(output_dir):
#         if base_name == '.DS_Store':
#             pass
#         elif not os.path.isdir(os.path.join(output_dir, base_name)):
#             pass
#         else:
#             base_dir = os.path.join(output_dir, base_name)
#             all_params = load_all_params(base_dir)
#             params = flatten_nested(all_params)
#             params.name = base_name

#             model_output_summary = model_output_summary.append(params)

#     return model_output_summary.T


# def filter_core(core, criteria, quality=None):
#     """
#     Filter simulations by criteria defined here
#     """
#     dum = core.copy()
#     for key, item in criteria.items():
#         dum = dum[dum[key] == item]
#     # if quality:
#     #     dum = dum[dum['quality'] > quality]
#     return dum
