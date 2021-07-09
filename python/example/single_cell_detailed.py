#!/usr/bin/env python3
# This script is included in documentation. Adapt line numbers if touched.

import arbor
import pandas
import seaborn
import sys
from arbor import mechanism as mech

# (1) Read the morphology from an SWC file.

# Read the SWC filename from input
# Example from docs: single_cell_detailed.swc

if len(sys.argv) < 2:
    print("No SWC file passed to the program")
    sys.exit(0)

filename = sys.argv[1]
morph = arbor.load_swc_arbor(filename)

# (2) Create and populate the label dictionary.

labels = arbor.label_dict()

# Regions:

# Add labels for tag 1, 2, 3, 4
labels['soma'] = '(tag 1)'
labels['axon'] = '(tag 2)'
labels['dend'] = '(tag 3)'
labels['last'] = '(tag 4)'
# Add a label for a region that includes the whole morphology
labels['all'] = '(all)'
# Add a label for the parts of the morphology with radius greater than 1.5 μm.
labels['gt_1.5'] = '(radius-ge (region "all") 1.5)'
# Join regions "last" and "gt_1.5"
labels['custom'] = '(join (region "last") (region "gt_1.5"))'

# Locsets:

# Add a labels for the root of the morphology and all the terminal points
labels['root']     = '(root)'
labels['terminal'] = '(terminal)'
# Add a label for the terminal locations in the "custom" region:
labels['custom_terminal'] = '(restrict (locset "terminal") (region "custom"))'
# Add a label for the terminal locations in the "axon" region:
labels['axon_terminal'] = '(restrict (locset "terminal") (region "axon"))'

# (3) Create and populate the decor.

decor = arbor.decor()

# Set the default properties of the cell (this overrides the model defaults).
decor.set_property(Vm =-55)
decor.set_ion('na', int_con=10,   ext_con=140, rev_pot=50, method='nernst/na')
decor.set_ion('k',  int_con=54.4, ext_con=2.5, rev_pot=-77)

# Override the cell defaults.
decor.paint('"custom"', tempK=270)
decor.paint('"soma"',   Vm=-50)

# Paint density mechanisms.
decor.paint('"all"', 'pas')
decor.paint('"custom"', 'hh')
decor.paint('"dend"',  mech('Ih', {'gbar': 0.001}))

# Place stimuli and spike detectors.
decor.place('"root"', arbor.iclamp(10, 1, current=2), 'iclamp0')
decor.place('"root"', arbor.iclamp(30, 1, current=2), 'iclamp1')
decor.place('"root"', arbor.iclamp(50, 1, current=2), 'iclamp2')
decor.place('"axon_terminal"', arbor.spike_detector(-10), 'detector')

# Single CV for the "soma" region
soma_policy = arbor.cv_policy_single('"soma"')
# Single CV for the "soma" region
dflt_policy = arbor.cv_policy_max_extent(1.0)
# default policy everywhere except the soma
policy = dflt_policy | soma_policy
# Set cv_policy
decor.discretization(policy)

# (4) Create the cell.

cell = arbor.cable_cell(morph, labels, decor)

# (5) Construct the model

model = arbor.single_cell_model(cell)

# (6) Set the model default properties

model.properties.set_property(Vm =-65, tempK=300, rL=35.4, cm=0.01)
model.properties.set_ion('na', int_con=10,   ext_con=140, rev_pot=50, method='nernst/na')
model.properties.set_ion('k',  int_con=54.4, ext_con=2.5, rev_pot=-77)

# Extend the default catalogue with the Allen catalogue.
# The function takes a second string parameter that can prefix
# the name of the mechanisms to avoid collisions between catalogues
# in this case we have no collisions so we use an empty prefix string.
model.catalogue.extend(arbor.allen_catalogue(), "")

# (7) Add probes.

# Add voltage probes on the "custom_terminal" locset
# which sample the voltage at 50 kHz
model.probe('voltage', where='"custom_terminal"',  frequency=50)

# (8) Run the simulation for 100 ms, with a dt of 0.025 ms

model.run(tfinal=100, dt=0.025)

# (9) Print the spikes.

print(len(model.spikes), 'spikes recorded:')
for s in model.spikes:
    print(s)

# (10) Plot the voltages

df = pandas.DataFrame()
for t in model.traces:
    df=df.append(pandas.DataFrame({'t/ms': t.time, 'U/mV': t.value, 'Location': str(t.location), 'Variable': t.variable}))

seaborn.relplot(data=df, kind="line", x="t/ms", y="U/mV",hue="Location",col="Variable",ci=None).savefig('single_cell_detailed_result.svg')
