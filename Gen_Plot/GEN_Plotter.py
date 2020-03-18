#!/usr/bin/env python
# 7 February 2019
# Abe Tishelman-Charny 
# Updated 31 May 2019 for HH MC 
# from ROOT import * 
import ROOT 
import array as arr 
import numpy as np
from GEN_Plotter_Config import * 
import subprocess
import ROOT
ROOT.gROOT.SetBatch(True)
from DataFormats.FWLite import Handle, Runs, Lumis, Events  #, ChainEvent 
import sys
import os 
from os import listdir
print '.'
print 'Plotting HH Variables'
print '.'
ol = '/uscms/home/rasharma/nobackup/double-higgs/HH_WWgg/Gen_Plot/myplots/'
# me, genHandle = -1, Handle('vector<reco::GenParticle>')
# summary = Handle('<LumiSummary>')
# For each variable
for iv,v in enumerate(vs):
    variable = v[0]
    print 'Plotting variable ', iv, ': ',v[0]
    # For each directory 
    histos = [] 
    # mg_title = v[0]
    # mg = ROOT.TMultiGraph('mg',mg_title)
    for dn,dir in enumerate(d):
        direc = dir[0] # full path of directory 
        rd = dir[1] # direc path with root://cmsxrootd.fnal.gov/ prefix for 'Events' module  
        mass = dir[2]
        channel = dir[3]  

        print "direc = ",direc
        print "rd = ",rd
        path_ends = [fp for fp in listdir(direc) if 'inLHE' not in fp]  # Save paths of non-LHE files
        #path_ends = [fp for fp in listdir(direc)]  # Save all paths 
        paths = []
        #print "path_ends = ",path_ends
        for pa in path_ends:
            tmp_path = rd + pa
            paths.append(tmp_path)
        # print '  Creating histos for file(s) of type: ',DID 
        # Get chosen particles to plot and number of them 
        pparams = []
        pparams = get_pparams(ptp)
        print "pparams = ",pparams
        xbins = v[1]
        xmin = v[2]
        xmax = v[3] # Mass_Channel_Variable 

        #Ha1 = ROOT.TLorentzVector()
        #Ha2 = ROOT.TLorentzVector()
        #Ha = ROOT.TLorentzVector()
        #HW1J1 = ROOT.TLorentzVector()
        #HW1J2 = ROOT.TLorentzVector()
        #HW1 = ROOT.TLorentzVector()
        #HW2J1 = ROOT.TLorentzVector()
        #HW2J2 = ROOT.TLorentzVector()
        #HW2 = ROOT.TLorentzVector()
        #HW = ROOT.TLorentzVector()



        plot_title = mass + ', ' + channel 
        h1 = ROOT.TH1F('h1',plot_title,xbins,xmin,xmax)
        h_pdgID_isHardProcess = ROOT.TH1F('h_pdgID_isHardProcess', 'particles from isHardProcess',30,0,30)
        h_pdgID_isHardProcess_status = ROOT.TH1F('h_pdgID_isHardProcess_status', 'particles from isHardProcess status',30,0,30)
        h_pdgID_isHardProcess_status_q = ROOT.TH1F('h_pdgID_isHardProcess_status_q', 'particles from isHardProcess status from q',30,0,30)
        h_pdgID_isHardProcess_status_g = ROOT.TH1F('h_pdgID_isHardProcess_status_g', 'particles from isHardProcess status from g',30,0,30)
        h_pdgID_isHardProcess_status_a = ROOT.TH1F('h_pdgID_isHardProcess_status_a', 'particles from isHardProcess status from a',30,0,30)
        h_pdgID_isHardProcess_status_w = ROOT.TH1F('h_pdgID_isHardProcess_status_w', 'particles from isHardProcess status from w',30,0,30)
        h_pdgID_isHardProcess_status_h = ROOT.TH1F('h_pdgID_isHardProcess_status_h', 'particles from isHardProcess status from h',30,0,30)

        h1_Photon0_pt = ROOT.TH1F("h1_Photon0_pt", ";Leading Photon p_{T} (GeV);",35,0,700)
        h1_Photon1_pt = ROOT.TH1F("h1_Photon1_pt", "",35,0,700)
        h1_Photon0_eta = ROOT.TH1F("h1_Photon0_eta", "",25,-5,5)
        h1_Photon1_eta = ROOT.TH1F("h1_Photon1_eta", "",25,-5,5)
        h1_Photon0_phi = ROOT.TH1F("h1_Photon0_phi", "",25,-5,5)
        h1_Photon1_phi = ROOT.TH1F("h1_Photon1_phi", "",25,-5,5)
        h1_Photon_dR = ROOT.TH1F("h1_Photon_dR", "",25,0,6)
        h2_Photon_dR_pT = ROOT.TH2F("h2_Photon_dR_pT", "", 25,0,6, 35,0,700)
        h1_Photon_InvMass = ROOT.TH1F("h1_Photon_InvMass", "",55,115,135)
        h1_WPlusJets0_pt = ROOT.TH1F("h1_WPlusJets0_pt", "",35,0,700)
        h1_WPlusJets1_pt = ROOT.TH1F("h1_WPlusJets1_pt", "",35,0,700)
        h1_WPlusJets0_eta = ROOT.TH1F("h1_WPlusJets0_eta", "",25,-5,5)
        h1_WPlusJets1_eta = ROOT.TH1F("h1_WPlusJets1_eta", "",25,-5,5)
        h1_WPlusJets0_phi = ROOT.TH1F("h1_WPlusJets0_phi", "",25,-5,5)
        h1_WPlusJets1_phi = ROOT.TH1F("h1_WPlusJets1_phi", "",25,-5,5)
        h1_WPlusJets_dR = ROOT.TH1F("h1_WPlusJets_dR", "",25,0,6)
        h2_WPlusJets_dR_pT = ROOT.TH2F("h2_WPlusJets_dR_pT", "", 25,0,6, 35,0,700)
        h1_WMinusJets0_pt = ROOT.TH1F("h1_WMinusJets0_pt", "",35,0,700)
        h1_WMinusJets1_pt = ROOT.TH1F("h1_WMinusJets1_pt", "",35,0,700)
        h1_WMinusJets0_eta = ROOT.TH1F("h1_WMinusJets0_eta", "",25,-5,5)
        h1_WMinusJets1_eta = ROOT.TH1F("h1_WMinusJets1_eta", "",25,-5,5)
        h1_WMinusJets0_phi = ROOT.TH1F("h1_WMinusJets0_phi", "",25,-5,5)
        h1_WMinusJets1_phi = ROOT.TH1F("h1_WMinusJets1_phi", "",25,-5,5)
        h1_WMinusJets_dR = ROOT.TH1F("h1_WMinusJets_dR", "",25,0,6)
        h2_WMinusJets_dR_pT = ROOT.TH2F("h2_WMinusJets_dR_pT", "", 25,0,6, 35,0,700)
        h1_WMinusCal_pt = ROOT.TH1F("h1_WMinusCal_pt", "",35,0,700)
        h1_WMinusCal_eta = ROOT.TH1F("h1_WMinusCal_eta", "",25,-5,5)
        h1_WMinusCal_phi = ROOT.TH1F("h1_WMinusCal_phi", "",25,-5,5)
        h1_WMinusCal_InvMass = ROOT.TH1F("h1_WMinusCal_InvMass", "",55,70,90)
        h1_WPlusCal_pt = ROOT.TH1F("h1_WPlusCal_pt", "",35,0,700)
        h1_WPlusCal_eta = ROOT.TH1F("h1_WPlusCal_eta", "",25,-5,5)
        h1_WPlusCal_phi = ROOT.TH1F("h1_WPlusCal_phi", "",25,-5,5)
        h1_WPlusCal_InvMass = ROOT.TH1F("h1_WPlusCal_InvMass", "",55,70,90)
        h1_WPlusWMinusCal_dR = ROOT.TH1F("h1_WPlusWMinusCal_dR", "",25,0,6)
        h2_WPlusWMinusCal_dR_pT = ROOT.TH2F("h2_WPlusWMinusCal_dR_pT", "", 25,0,6, 35,0,700)
        h1_WMinus_pt = ROOT.TH1F("h1_WMinus_pt", "",35,0,700)
        h1_WMinus_eta = ROOT.TH1F("h1_WMinus_eta", "",25,-5,5)
        h1_WMinus_phi = ROOT.TH1F("h1_WMinus_phi", "",25,-5,5)
        h1_WMinus_InvMass = ROOT.TH1F("h1_WMinus_InvMass", "",55,70,90)
        h1_WPlus_pt = ROOT.TH1F("h1_WPlus_pt", "",35,0,700)
        h1_WPlus_eta = ROOT.TH1F("h1_WPlus_eta", "",25,-5,5)
        h1_WPlus_phi = ROOT.TH1F("h1_WPlus_phi", "",25,-5,5)
        h1_WPlus_InvMass = ROOT.TH1F("h1_WPlus_InvMass", "",55,70,90)
        h1_WPlusWminus_dR = ROOT.TH1F("h1_WPlusWminus_dR", "",25,0,6)
        h2_WPlusWminus_dR_pT = ROOT.TH2F("h2_WPlusWminus_dR_pT", "", 25,0,6, 35,0,700)
        h1_WPlusWminus_M = ROOT.TH1F("h1_WPlusWminus_M", "",55,115,135)

        count_events = 0
        for ip,path in enumerate(paths):
            if ip == max_files: break 
            print '    Processing file ', ip+1, ': ',path 
            events = Events(path) # needs to be file with root prefix                 
            # db = (float(xmax) - float(xmin)) / float(xbins) 
            # Loop events 
            # Add to sum 
            for iev, event in enumerate(events):
                count_events += 1
                if iev%500 == 0: print'on event',iev 
                if iev == me: break # Max events 
                #print "event number = ",iev
                
                #event.getByLabel('prunedGenParticles', genHandle)
                events.getByLabel('genParticles', genHandle)
                genParticles = genHandle.product()
                events.getByLabel('ak4GenJets', ak4genHandle)
                ak4GenJets = ak4genHandle.product()
                events.getByLabel('ak8GenJets', ak8genHandle)
                ak8GenJets = ak4genHandle.product()

                # Fill histograms with current variable 
                for params in pparams:
                    particle = params[0]
                    nparticles = params[1]
                    pdgIDs = params[2]
                    #print "particle = ",particle,"\tnparticles = ",nparticles,"\tpdgIDs = ",pdgIDs
                    #ps = [p for p in genParticles if p.isHardProcess() and abs(p.pdgId()) in pdgIDs]
                    #print "="*30
                    
                    PhotonFromHiggs = []
                    QuarkFromWp = []
                    QuarkFromWm = []
                    WFromHiggs = []
                    Higgs = []

                    for p in genParticles:
                       if p.isHardProcess():
                          h_pdgID_isHardProcess.Fill(p.pdgId())
                          h_pdgID_isHardProcess_status.Fill(p.status())
                          #print type(p.p4())
                          if abs(p.pdgId()) < 5: 
                             h_pdgID_isHardProcess_status_q.Fill(p.status())
                          if abs(p.pdgId()) == 21: 
                             h_pdgID_isHardProcess_status_g.Fill(p.status())
                          if abs(p.pdgId()) == 22: 
                             h_pdgID_isHardProcess_status_a.Fill(p.status())
                          if abs(p.pdgId()) == 24: 
                             h_pdgID_isHardProcess_status_w.Fill(p.status())
                          if abs(p.pdgId()) == 25: 
                             h_pdgID_isHardProcess_status_h.Fill(p.status())
                          #if p.pdgId() == -24: print p.pdgId(),'\t',p.daughter(0).pdgId(),'\t',p.status(),'\t',p.daughter(0).status()
                          #if p.pdgId() == 24: print p.pdgId(),'\t',p.daughter(0).pdgId(),'\t',p.status(),'\t',p.daughter(0).status()

                          if ( abs(p.pdgId()) == 22 and 
                               p.mother(0).pdgId() == 25 and 
                               (p.status() == 23 or p.status() == 1)
                             ): 
                             PhotonFromHiggs.append(p.p4())
                          if ( abs(p.pdgId()) < 5 and 
                               p.mother(0).pdgId() == 24 and 
                               p.status() == 23
                             ): 
                             QuarkFromWp.append(p.p4())
                          if ( abs(p.pdgId()) < 5 and 
                               p.mother(0).pdgId() == -24 and 
                               p.status() == 23
                             ): 
                             QuarkFromWm.append(p.p4())
                          if ( abs(p.pdgId()) == 24 and 
                               abs(p.mother(0).pdgId()) == 25 and 
                               p.status() == 22 and
                               abs(p.daughter(0).pdgId()) < 5
                             ): 
                             WFromHiggs.append(p.p4())
                          if ( abs(p.pdgId()) == 25 and
                               #( abs(p.daughter(0).pdgId()) == 22 or
                               #  abs(p.daughter(0).pdgId()) == 24
                               #) and
                               p.status() == 22
                             ):
                             Higgs.append(p.p4())
                             #print "daughter = ",abs(p.daughter(0).pdgId()),"\tStatus = ",p.daughter(0).status()
                          """
                          if ( abs(p.pdgId()) == 25 and 
                               p.status() == 22 and 
                               ( ( p.daughter(0).pdgId() == 22 and 
                                   ( p.daughter(0).status() == 23 or p.daughter(0).status() == 1 ) 
                                 ) or 
                                 ( abs(p.daughter(0).pdgId())==24 and 
                                   abs(p.daughter(0).daughter(0).pdgId()) < 5  
                                 )
                               )
                              ): 
                              Higgs.append(p.p4())
                          """
                    #print "*"*30
                    #print "Size: Photons = ",len(PhotonFromHiggs)
                    #print "Size: QuarkFromWp = ",len(QuarkFromWp)
                    #print "Size: QuarkFromWm = ",len(QuarkFromWm)
                    #print "Size: WFromHiggs = ",len(WFromHiggs)
                    #print "Size: Higgs = ",len(Higgs)
                    ##print "Higgs pt = ",Higgs[0].pt(),'\t',Higgs[1].pt()
                    ##print "Photons pt = ",PhotonFromHiggs[0].pt(),'\t',PhotonFromHiggs[1].pt()
                    if not (len(PhotonFromHiggs) == 2 and len(QuarkFromWp) == 2 and len(QuarkFromWm) == 2 and len(WFromHiggs) == 2 and len(Higgs) == 2): 
                      continue

                    Higgs = order_particle_4vec(Higgs)
                    PhotonFromHiggs = order_particle_4vec(PhotonFromHiggs)
                    QuarkFromWp = order_particle_4vec(QuarkFromWp)
                    QuarkFromWm = order_particle_4vec(QuarkFromWm)
                    WFromHiggs = order_particle_4vec(WFromHiggs)
                    
                    #Fill photon variables
                    h1_Photon0_pt.Fill(PhotonFromHiggs[0].pt())
                    h1_Photon1_pt.Fill(PhotonFromHiggs[1].pt())
                    h1_Photon0_eta.Fill(PhotonFromHiggs[0].eta())
                    h1_Photon1_eta.Fill(PhotonFromHiggs[1].eta())
                    h1_Photon0_phi.Fill(PhotonFromHiggs[0].phi())
                    h1_Photon1_phi.Fill(PhotonFromHiggs[1].phi())
                    h1_Photon_dR.Fill( ROOT.Math.VectorUtil.DeltaR(PhotonFromHiggs[0], PhotonFromHiggs[1]) )
                    h2_Photon_dR_pT.Fill(ROOT.Math.VectorUtil.DeltaR(PhotonFromHiggs[0], PhotonFromHiggs[1]), (PhotonFromHiggs[0]+PhotonFromHiggs[1]).pt()) 
                    h1_Photon_InvMass.Fill( (PhotonFromHiggs[0] + PhotonFromHiggs[1]).M())


                    #Fill jets from w+ variables
                    h1_WPlusJets0_pt.Fill(QuarkFromWp[0].pt())
                    h1_WPlusJets1_pt.Fill(QuarkFromWp[1].pt())
                    h1_WPlusJets0_eta.Fill(QuarkFromWp[0].eta())
                    h1_WPlusJets1_eta.Fill(QuarkFromWp[1].eta())
                    h1_WPlusJets0_phi.Fill(QuarkFromWp[0].phi())
                    h1_WPlusJets1_phi.Fill(QuarkFromWp[1].phi())
                    h1_WPlusJets_dR.Fill( ROOT.Math.VectorUtil.DeltaR(QuarkFromWp[0], QuarkFromWp[1]) )
                    h2_WPlusJets_dR_pT.Fill( ROOT.Math.VectorUtil.DeltaR(QuarkFromWp[0], QuarkFromWp[1]), (QuarkFromWp[0] + QuarkFromWp[1]).pt())

                    #Fill jets from w- variables
                    h1_WMinusJets0_pt.Fill(QuarkFromWm[0].pt())
                    h1_WMinusJets1_pt.Fill(QuarkFromWm[1].pt())
                    h1_WMinusJets0_eta.Fill(QuarkFromWm[0].eta())
                    h1_WMinusJets1_eta.Fill(QuarkFromWm[1].eta())
                    h1_WMinusJets0_phi.Fill(QuarkFromWm[0].phi())
                    h1_WMinusJets1_phi.Fill(QuarkFromWm[1].phi())
                    h1_WMinusJets_dR.Fill( ROOT.Math.VectorUtil.DeltaR(QuarkFromWm[0], QuarkFromWm[1]) )
                    h2_WMinusJets_dR_pT.Fill(ROOT.Math.VectorUtil.DeltaR(QuarkFromWm[0], QuarkFromWm[1]), (QuarkFromWm[0]+QuarkFromWm[1]).pt())

                    # Fill W's attribute calculated from 4-vector
                    h1_WMinusCal_pt.Fill( (QuarkFromWp[0]+QuarkFromWp[1]).pt() )
                    h1_WMinusCal_eta.Fill( (QuarkFromWp[0]+QuarkFromWp[1]).eta() )
                    h1_WMinusCal_phi.Fill( (QuarkFromWp[0]+QuarkFromWp[1]).phi() )
                    h1_WMinusCal_InvMass.Fill( (QuarkFromWp[0]+QuarkFromWp[1]).M() )

                    h1_WPlusCal_pt.Fill( (QuarkFromWm[0]+QuarkFromWm[1]).pt() )
                    h1_WPlusCal_eta.Fill( (QuarkFromWm[0]+QuarkFromWm[1]).eta() )
                    h1_WPlusCal_phi.Fill( (QuarkFromWm[0]+QuarkFromWm[1]).phi() )
                    h1_WPlusCal_InvMass.Fill((QuarkFromWm[0]+QuarkFromWm[1]).M() )

                    h1_WPlusWMinusCal_dR.Fill( ROOT.Math.VectorUtil.DeltaR((QuarkFromWm[0]+QuarkFromWm[1]), (QuarkFromWp[0]+QuarkFromWp[1])) )
                    h2_WPlusWMinusCal_dR_pT.Fill(ROOT.Math.VectorUtil.DeltaR((QuarkFromWm[0]+QuarkFromWm[1]), (QuarkFromWp[0]+QuarkFromWp[1])), (QuarkFromWm[0]+QuarkFromWm[1]+QuarkFromWp[0]+QuarkFromWp[1]).pt())
                    
                    # Fill W's attribute access directly
                    h1_WMinus_pt.Fill(WFromHiggs[0].pt())
                    h1_WMinus_eta.Fill(WFromHiggs[0].eta())
                    h1_WMinus_phi.Fill(WFromHiggs[0].phi())
                    h1_WMinus_InvMass.Fill(WFromHiggs[0].M())

                    h1_WPlus_pt.Fill(WFromHiggs[1].pt())
                    h1_WPlus_eta.Fill(WFromHiggs[1].eta())
                    h1_WPlus_phi.Fill(WFromHiggs[1].phi())
                    h1_WPlus_InvMass.Fill(WFromHiggs[1].M())
                    h1_WPlusWminus_dR.Fill( ROOT.Math.VectorUtil.DeltaR(WFromHiggs[0], WFromHiggs[1]) )
                    h2_WPlusWminus_dR_pT.Fill(ROOT.Math.VectorUtil.DeltaR(WFromHiggs[0], WFromHiggs[1]), (WFromHiggs[0]+WFromHiggs[1]).pt())
                    h1_WPlusWminus_M.Fill( (WFromHiggs[0]+WFromHiggs[1]).M())
                    

                    

                    #        val = invmass(ps[0].p4(),ps[1].p4())
        output_path = ol + mass + '_' + channel + '_' + particle + '_' + variable 
        c1 = ROOT.TCanvas()
        color = dir[4]
        print 'h_pdgID_isHardProcess entry = ',h_pdgID_isHardProcess.GetEntries()

        h_pdgID_isHardProcess.Draw();  c1.SaveAs("h_pdgID_isHardProcess.png")
        h_pdgID_isHardProcess_status.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status.png")
        h_pdgID_isHardProcess_status_q.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status_q.png")
        h_pdgID_isHardProcess_status_g.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status_g.png")
        h_pdgID_isHardProcess_status_a.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status_a.png")
        h_pdgID_isHardProcess_status_w.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status_w.png")
        h_pdgID_isHardProcess_status_h.Draw();  c1.SaveAs("h_pdgID_isHardProcess_status_h.png")
        h1_Photon0_pt.Draw();  c1.SaveAs("h1_Photon0_pt.png")
        h1_Photon1_pt.Draw();  c1.SaveAs("h1_Photon1_pt.png")
        h1_Photon0_eta.Draw();  c1.SaveAs("h1_Photon0_eta.png")
        h1_Photon1_eta.Draw();  c1.SaveAs("h1_Photon1_eta.png")
        h1_Photon0_phi.Draw();  c1.SaveAs("h1_Photon0_phi.png")
        h1_Photon1_phi.Draw();  c1.SaveAs("h1_Photon1_phi.png")
        h1_Photon_dR.Draw();  c1.SaveAs("h1_Photon_dR.png")
        h2_Photon_dR_pT.Draw("colz"); c1.SaveAs("h2_Photon_dR_pT.png")
        h1_Photon_InvMass.Draw();  c1.SaveAs("h1_Photon_InvMass.png")
        h1_WPlusJets0_pt.Draw();  c1.SaveAs("h1_WPlusJets0_pt.png")
        h1_WPlusJets1_pt.Draw();  c1.SaveAs("h1_WPlusJets1_pt.png")
        h1_WPlusJets0_eta.Draw();  c1.SaveAs("h1_WPlusJets0_eta.png")
        h1_WPlusJets1_eta.Draw();  c1.SaveAs("h1_WPlusJets1_eta.png")
        h1_WPlusJets0_phi.Draw();  c1.SaveAs("h1_WPlusJets0_phi.png")
        h1_WPlusJets1_phi.Draw();  c1.SaveAs("h1_WPlusJets1_phi.png")
        h1_WPlusJets_dR.Draw();  c1.SaveAs("h1_WPlusJets_dR.png")
        h2_WPlusJets_dR_pT.Draw("colz"); c1.SaveAs("h2_WPlusJets_dR_pT.png")
        h1_WMinusJets0_pt.Draw();  c1.SaveAs("h1_WMinusJets0_pt.png")
        h1_WMinusJets1_pt.Draw();  c1.SaveAs("h1_WMinusJets1_pt.png")
        h1_WMinusJets0_eta.Draw();  c1.SaveAs("h1_WMinusJets0_eta.png")
        h1_WMinusJets1_eta.Draw();  c1.SaveAs("h1_WMinusJets1_eta.png")
        h1_WMinusJets0_phi.Draw();  c1.SaveAs("h1_WMinusJets0_phi.png")
        h1_WMinusJets1_phi.Draw();  c1.SaveAs("h1_WMinusJets1_phi.png")
        h1_WMinusJets_dR.Draw();  c1.SaveAs("h1_WMinusJets_dR.png")
        h2_WMinusJets_dR_pT.Draw("colz"); c1.SaveAs("h2_WMinusJets_dR_pT.png")
        h1_WMinusCal_pt.Draw();  c1.SaveAs("h1_WMinusCal_pt.png")
        h1_WMinusCal_eta.Draw();  c1.SaveAs("h1_WMinusCal_eta.png")
        h1_WMinusCal_phi.Draw();  c1.SaveAs("h1_WMinusCal_phi.png")
        h1_WMinusCal_InvMass.Draw();  c1.SaveAs("h1_WMinusCal_InvMass.png")
        h1_WPlusCal_pt.Draw();  c1.SaveAs("h1_WPlusCal_pt.png")
        h1_WPlusCal_eta.Draw();  c1.SaveAs("h1_WPlusCal_eta.png")
        h1_WPlusCal_phi.Draw();  c1.SaveAs("h1_WPlusCal_phi.png")
        h1_WPlusCal_InvMass.Draw();  c1.SaveAs("h1_WPlusCal_InvMass.png")
        h1_WPlusWMinusCal_dR.Draw();  c1.SaveAs("h1_WPlusWMinusCal_dR.png")
        h2_WPlusWMinusCal_dR_pT.Draw("colz"); c1.SaveAs("h2_WPlusWMinusCal_dR_pT.png")
        h1_WMinus_pt.Draw();  c1.SaveAs("h1_WMinus_pt.png")
        h1_WMinus_eta.Draw();  c1.SaveAs("h1_WMinus_eta.png")
        h1_WMinus_phi.Draw();  c1.SaveAs("h1_WMinus_phi.png")
        h1_WMinus_InvMass.Draw();  c1.SaveAs("h1_WMinus_InvMass.png")
        h1_WPlus_pt.Draw();  c1.SaveAs("h1_WPlus_pt.png")
        h1_WPlus_eta.Draw();  c1.SaveAs("h1_WPlus_eta.png")
        h1_WPlus_phi.Draw();  c1.SaveAs("h1_WPlus_phi.png")
        h1_WPlus_InvMass.Draw();  c1.SaveAs("h1_WPlus_InvMass.png")
        h1_WPlusWminus_dR.Draw();  c1.SaveAs("h1_WPlusWminus_dR.png")
        h2_WPlusWminus_dR_pT.Draw("colz"); c1.SaveAs("h2_WPlusWminus_dR_pT.png")
        h1_WPlusWminus_M.Draw(); c1.SaveAs("h1_WPlusWminus_M.png")

        histos.append(h_pdgID_isHardProcess)
        print "total events = ",count_events

    print'histos = ',histos 

    # c2 = ROOT.TCanvas()
    # for ih,h in enumerate(histos):
    #     print'h = ',h
    #     h.SetStats(0)
    #     if ih == 0:
    #         h.Draw()
    #     else:
    #         h.Draw('same')
    # xmin, ymin, xmax, ymax = 0.6,0.5,0.8,0.7
    # legend = ROOT.TLegend(xmin,ymin,xmax,ymax)
    # for ih,h in enumerate(histos):
    #     legend.AddEntry(h)
    # legend.Draw('same')
    # c2.Update()
    # c2.SaveAs(ol + variable + '_Combined.png')
                                   
