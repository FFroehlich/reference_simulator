import os
import sys
import importlib
import amici
import pandas as pd
import numpy as np

model_name = 'ERBB_RAS_AKT_Drugs'

sys.path.insert(0, os.path.abspath(model_name))
model_module = importlib.import_module(model_name)

model = model_module.getModel()
# run model to steadystate
model.setTimepoints([np.infty])

solver = model.getSolver()

conditions = pd.read_csv('conditions_petab.tsv', sep='\t')


def run_simulation(
        drug_concs,
        cell_line='A2058',
):
    condition = conditions.loc[
        conditions.conditionId == f'TUMOR-{cell_line}-cellline-01-01', :
    ]

    if len(condition) == 0:
        raise ValueError(
            f'Requested cell-line "{cell_line}" has no condition data.'
        )

    # set the fixed parameters
    for col in condition.columns:
        if col in model.getFixedParameterIds():
            model.setFixedParameterById(col, condition[col].values[0])

    edata_ref = amici.ExpData(model.get())

    for drug, conc in drug_concs.items():
        parameter_id = model.getFixedParameterIds()[
            model.getFixedParameterNames().index(drug)
        ]
        model.setFixedParameterById(parameter_id, conc)
        print(f'{drug}: {conc}')

    # extract simulation conditions that were specified in the model
    edata_cond = amici.ExpData(model.get())
    edatas = [edata_ref, edata_cond]
    rdatas = amici.runAmiciSimulations(model, solver, edatas)
    print(f'time to steadystate {rdatas[0]["t_steadystate"]}')
    print(f'relative proliferation '
          f'{rdatas[1]["y"][0, 0]/rdatas[0]["y"][0, 0]}')
    return amici.getSimulationStatesAsDataFrame(model, edatas, rdatas)


concentrations = {
    'PD0325901': 0.0,
    'PLX-4720': 0.0,
    'Selumetinib': 0.0,
    'Lapatinib': 0.0,
    'Erlotinib': 0.0,
    'CHIR-265': 0.0,
    'Vandetanib': 0.0,
}

df = run_simulation(drug_concs=concentrations,
                    cell_line='A2058')

df.to_csv('simulation_results.csv')

