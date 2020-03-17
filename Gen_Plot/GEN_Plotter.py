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
        plot_title = mass + ', ' + channel 
        h1 = ROOT.TH1F('h1',plot_title,xbins,xmin,xmax)
        h_pdgID_isHardProcess = ROOT.TH1F('h_pdgID_isHardProcess', 'particles from isHardProcess',30,0,30)
        h_pdgID_isHardProcess_status = ROOT.TH1F('h_pdgID_isHardProcess_status', 'particles from isHardProcess',30,0,30)
        h_pdgID_isHardProcess_status_q = ROOT.TH1F('h_pdgID_isHardProcess_status_q', 'particles from isHardProcess from q',30,0,30)
        h_pdgID_isHardProcess_status_g = ROOT.TH1F('h_pdgID_isHardProcess_status_g', 'particles from isHardProcess from g',30,0,30)
        h_pdgID_isHardProcess_status_a = ROOT.TH1F('h_pdgID_isHardProcess_status_a', 'particles from isHardProcess from a',30,0,30)
        h_pdgID_isHardProcess_status_w = ROOT.TH1F('h_pdgID_isHardProcess_status_w', 'particles from isHardProcess from w',30,0,30)
        h_pdgID_isHardProcess_status_h = ROOT.TH1F('h_pdgID_isHardProcess_status_h', 'particles from isHardProcess from h',30,0,30)
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
                    print "="*30
                    for p in genParticles:
                       if p.isHardProcess():
                          h_pdgID_isHardProcess.Fill(p.pdgId())
                          h_pdgID_isHardProcess_status.Fill(p.status())
                          if abs(p.pdgId()) < 5: h_pdgID_isHardProcess_status_q.Fill(p.status())
                          if abs(p.pdgId()) == 21: h_pdgID_isHardProcess_status_g.Fill(p.status())
                          if abs(p.pdgId()) == 22: h_pdgID_isHardProcess_status_a.Fill(p.status())
                          if abs(p.pdgId()) == 24: h_pdgID_isHardProcess_status_w.Fill(p.status())
                          if abs(p.pdgId()) == 25: h_pdgID_isHardProcess_status_h.Fill(p.status())
                          if p.pdgId() == -24: print p.daughter(0).pdgId(),"\t",p.daughter(1).pdgId()
                          if p.pdgId() == 24: print p.daughter(0).pdgId(),"\t",p.daughter(1).pdgId()
                          if abs(p.pdgId()) < 5:
                             print p.pdgId(),"\t",p.mother(0).pdgId()
                          #print p.pdgId()
                       #if p.isHardProcess() and abs(p.pdgId()) == 24: 
                       #if abs(p.pdgId()) == 24: 
                       #   print p.pdgId(),"\t",p.status(),"\t",type(p)
                    #print "="*30
                    #for p in ak4GenJets:
                    #   print p.pt()
                    #print "="*30
                    #for p in ak8GenJets:
                    #   print p.pt()
                    #ps = [p for p in genParticles if p.isHardProcess() and abs(p.pdgId()) in pdgIDs and abs(p.daughter(0).pdgId() == 22)]   
                    #print "particle = ",ps[0].pdgId(),"\t",ps[1].pdgId()
                    #print "particles = ",type(ps),"\t",ps[0].pdgId()
                    #for i, pp in enumerate(ps):
                    #   print type(pp),"\t",pp.pdgId()
                    #for i,p in enumerate(ps): print p[i].pdgId(),"\t",
                    #print "\n\n.....\n"

                    # ps = [p for p in genParticles if p.isHardProcess() and abs(p.pdgId()) in pdgIDs and abs(p.daughter(0).pdgId() == 25)]   
                    # ps = [p for p in genParticles if p.isHardProcess() and abs(p.pdgId()) in pdgIDs]   
                    # for p in ps:
                    #     print'p = ',p.p4()
                    #if v[0] == 'pt':
                    #    val = ps[0].p4().pt() 

                    #if nparticles == 4:
                    #    if v[0] == 'invm':
                    #        val = invmass(ps[0].p4(),ps[1].p4())
                    #        # val = ps[0].p4().pt()
                    #        # if particle == 'R':
                    #            # avoid double count 
                    #        h1.Fill(val)
                    #        # get invmass 
                    #else: 
                    #    for p in ps:
                    #        #val = eval("p." + v[0] + "()")
                    #        val = ps[0].p4().pt()
                    #        #h1.Fill(val)
        output_path = ol + mass + '_' + channel + '_' + particle + '_' + variable 
        c1 = ROOT.TCanvas()
        color = dir[4]
        h_pdgID_isHardProcess.SetLineColor(eval(color))
        h_pdgID_isHardProcess.Draw() 
        c1.SaveAs("h_pdgID_isHardProcess" + ".png")
        h_pdgID_isHardProcess.SaveAs("h_pdgID_isHardProcess" + ".C")
        h_pdgID_isHardProcess.SaveAs("h_pdgID_isHardProcess" + ".root")    
        c1.SetLogy(1)
        c1.SaveAs("h_pdgID_isHardProcess" + "Log.png")
        c1.SetLogy(0)
        h_pdgID_isHardProcess_status.SetLineColor(eval(color))
        h_pdgID_isHardProcess_status.Draw() 
        c1.SaveAs("h_pdgID_isHardProcess_status" + ".png")
        h_pdgID_isHardProcess_status.SaveAs("h_pdgID_isHardProcess_status" + ".C")
        h_pdgID_isHardProcess_status.SaveAs("h_pdgID_isHardProcess_status" + ".root")    
        c1.SetLogy(1)
        c1.SaveAs("h_pdgID_isHardProcess_status" + "Log.png")
        c1.SetLogy(0)
        h_pdgID_isHardProcess_status_q.Draw()
        c1.SaveAs("h_pdgID_isHardProcess_status_q.png")
        h_pdgID_isHardProcess_status_g.Draw()
        c1.SaveAs("h_pdgID_isHardProcess_status_g.png")
        h_pdgID_isHardProcess_status_a.Draw()
        c1.SaveAs("h_pdgID_isHardProcess_status_a.png")
        h_pdgID_isHardProcess_status_w.Draw()
        c1.SaveAs("h_pdgID_isHardProcess_status_w.png")
        h_pdgID_isHardProcess_status_h.Draw()
        c1.SaveAs("h_pdgID_isHardProcess_status_h.png")
        print 'h_pdgID_isHardProcess entry = ',h_pdgID_isHardProcess.GetEntries()
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
                                   
