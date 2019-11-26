from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import pyjet
import fastjet as fj
from pyjet import cluster, DTYPE_PTEPM
from pyjet.testdata import get_event
from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap

qcd_data = pd.read_hdf('events_LHCO2020_backgroundMC_Pythia.h5')
sig_data = pd.read_hdf('events_LHCO2020_BlackBox1.h5')

def run(data, n_events = 1000000):
    
    out = []
    njets = []
    met = []
    rho = []
    jets_corr = []
    area = []

    # run jet clustering with AntiKt, R=1.0
    R = 1.0
    jet_def = fj.JetDefinition(fj.antikt_algorithm, R)

    # Area Definition
    ghost_maxrap = 4.7
    area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(ghost_maxrap))

    # Background estimator
    select_rapidity = fj.SelectorAbsRapMax(ghost_maxrap)
    bge = fj.JetMedianBackgroundEstimator(select_rapidity, jet_def, area_def)

    # Loop over events
    for ievt in range(n_events):

        # Build a list of all particles and MET Information
        pjs = []
        pjmet = fj.PseudoJet()
        pjmets = fj.PseudoJet()

        for i in range(int(data.shape[1]/3)):

            pj = fj.PseudoJet()
            pj.reset_PtYPhiM(data.at[ievt, 3*i+0], data.at[ievt, 3*i+1], data.at[ievt, 3*i+2], 0)
            if pj.pt() > 1.0:
                pjs.append(pj)

            pjmet.reset_momentum(-pj.px(), -pj.py(), 0, pj.E())
            pjmets = pjmets + pjmet

        met.append(pjmets)

        # Cluster sequence
        clust_seq = fj.ClusterSequenceArea(pjs, jet_def, area_def)
        jets = fj.sorted_by_pt(clust_seq.inclusive_jets())
        jets = [j for j in jets if j.pt() > 30. and j.eta() < 4.7]
        area.append([jets[0].area(), jets[1].area()])

        # Save the two leading jets and njets
        jets_sel = (fj.SelectorNHardest(2))(jets)
        out.append(jets_sel)
        njets.append(len(jets))

        # Background estimator
        bge.set_particles(pjs)
        rho.append(bge.rho())

        # Correct Jets
        sub = fj.Subtractor(bge)
        sub_jets = sub(jets_sel)
        jets_corr.append(sub_jets)

    return out, area, met, njets, rho, jets_corr

# Outputs 
out_qcd, area_qcd, met_qcd, njets_qcd, rho_qcd, corr_qcd = run(qcd_data)
out_sig, area_sig, met_sig, njets_sig, rho_sig, corr_sig = run(sig_data)

# Histos Initialization
key_names = ['mjj', 'j1pt', 'j1eta', 'j1phi', 'j2pt', 'j2eta', 'j2phi', 'j1area', 'j2area', 'met_pt', 'met_phi', 'njets', 'rho', 'mjj_corr', 'j1pt_corr', 'j2pt_corr']
histos_sig = {}
histos_qcd = {}
for k in key_names:
    histos_sig[k] = []
    histos_qcd[k] = []

# Histos for Signal
for evt in out_sig:
    histos_sig['mjj'].append((evt[0]+evt[1]).m())
    histos_sig['j1pt'].append(evt[0].pt())
    histos_sig['j1eta'].append(evt[0].eta())
    histos_sig['j1phi'].append(evt[0].phi())
    histos_sig['j2pt'].append(evt[1].pt())
    histos_sig['j2eta'].append(evt[1].eta())
    histos_sig['j2phi'].append(evt[1].phi())
for evt in area_sig:
    histos_sig['j1area'].append(evt[0])
    histos_sig['j2area'].append(evt[1])
for evt in met_sig:
    histos_sig['met_pt'].append(evt.pt())
    histos_sig['met_phi'].append(evt.phi())
for evt in njets_sig:
    histos_sig['njets'].append(evt)
for evt in rho_sig:
    histos_sig['rho'].append(evt)
for evt in corr_sig:
    histos_sig['mjj_corr'].append((evt[0]+evt[1]).m())
    histos_sig['j1pt_corr'].append(evt[0].pt())
    histos_sig['j2pt_corr'].append(evt[1].pt())

# Histos for QCD
for evt in out_qcd:
    histos_qcd['mjj'].append((evt[0]+evt[1]).m())
    histos_qcd['j1pt'].append(evt[0].pt())
    histos_qcd['j1eta'].append(evt[0].eta())
    histos_qcd['j1phi'].append(evt[0].phi())
    histos_qcd['j2pt'].append(evt[1].pt())
    histos_qcd['j2eta'].append(evt[1].eta())
    histos_qcd['j2phi'].append(evt[1].phi())
for evt in area_qcd:
    histos_qcd['j1area'].append(evt[0])
    histos_qcd['j2area'].append(evt[1])
for evt in met_qcd:
    histos_qcd['met_pt'].append(evt.pt())
    histos_qcd['met_phi'].append(evt.phi())
for evt in njets_qcd:
    histos_qcd['njets'].append(evt)
for evt in rho_qcd:
    histos_qcd['rho'].append(evt)
for evt in corr_qcd:
    histos_qcd['mjj_corr'].append((evt[0]+evt[1]).m())
    histos_qcd['j1pt_corr'].append(evt[0].pt())
    histos_qcd['j2pt_corr'].append(evt[1].pt())

# Plot the histos and save
xlabel_names = ['Di-jet Invariant Mass (GeV)', 'Jet 1 pT (GeV)', 'Jet 1 Eta', 'Jet 1 Phi', 'Jet 2 pT (GeV)', 'Jet 2 Eta', 'Jet 2 Phi', 'Jet 1 Area', 'Jet 2 Area', 'MET pT (GeV)', 'MET Phi', 'Number of Jets', 'Energy Density', 'Corrected Di-jet Invariant Mass (GeV)', 'Corrected Jet 1 pT (GeV)', 'Corrected Jet 2 pT (GeV)']

for i, pn in enumerate(key_names):
    _ = plt.hist(histos_sig[pn], label="Signal", bins=60, histtype='step')
    _ = plt.hist(histos_qcd[pn], label="QCD", bins=60, histtype='step')
    plt.legend()
    plt.title('QCD vs Signal')
    plt.xlabel(xlabel_names[i])
    plt.ylabel('Events per bin')
    plt.savefig('plot/'+pn+'.png')
    plt.clf()
