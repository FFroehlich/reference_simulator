import libsbml
import importlib
import amici
import os
import sys
import pandas as pd
import petab.sbml

# SBML model we want to import
sbml_file = 'CS_Signalling_ERBB_RAS_AKT_petab.xml'
# Name of the model that will also be the name of the python module
model_name = 'ERBB_RAS_AKT_Drugs'
# Directory to which the generated model code is written
model_output_dir = model_name

sbml_importer = amici.SbmlImporter(sbml_file)

# extract observable definition from sbml
observables = amici.assignmentRules2observables(
    sbml_importer.sbml,
    filter_function=lambda p: p.getName() in ['observable_proliferation']
)

petab.sbml.constant_species_to_parameters(sbml_importer.sbml)

libsbml.writeSBMLToFile(sbml_importer.sbml_doc,
                        'CS_Signalling_ERBB_RAS_AKT_modified.xml')

condition_table = pd.read_csv('conditions_petab.tsv', sep='\t')

# condition parameters should be everything that is defined in conditions and
# also specified in the model
constantParameters = [
    par for par in condition_table.columns
    if sbml_importer.sbml.getParameter(par)
]

sbml_importer.sbml2amici(model_name,
                         model_output_dir,
                         verbose=False,
                         observables=observables,
                         constantParameters=constantParameters)

sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)

model = model_module.getModel()

print("Model parameters:", model.getParameterIds())
print("Model outputs:   ", model.getObservableIds())
print("Model states:    ", model.getStateIds())
