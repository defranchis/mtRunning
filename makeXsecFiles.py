
import os
import rundec
import numpy as np

import ROOT as rt
import variables as var
import constants as cnst

from ROOT import TString, TH2D, TRandom3, TF1, TGraph, TLine, TCanvas, TGraphErrors, TLegend, TLatex
from variables import xsec_1, xsec_2, xsec_3


# options

estimate_contribs = True
ntoys = 100

replace_corr = False

# end options



################################

# "setMasses" initializes which masses should be considered

################################

def setMasses():

    global mass_v
    mass_v = []
    i_mass = cnst.mass_max
    while i_mass >= cnst.mass_min :
        mass_v.append(i_mass)
        if i_mass <= cnst.mass_fine_max and i_mass > cnst.mass_fine_min+.1 : i_mass-= 0.2
        else : i_mass-= 0.5
    return


################################

# "mtmt2mtmu" converts mt(mt) to mt(mu) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmt2mtmu(mt, mu):

    crd = rundec.CRunDec()

    asmt = crd.AlphasExact(cnst.asMZ, cnst.MZ, mt, cnst.nflav, cnst.nloops)
    asmu = crd.AlphasExact(cnst.asMZ, cnst.MZ, mu, cnst.nflav, cnst.nloops)
    
    mtmu = crd.mMS2mMS(mt, asmt, asmu, cnst.nflav, cnst.nloops)
    
    return mtmu



################################

# "mtmu2mtmu" converts mt(mu1) to mt(mu2) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmu2mtmu (mtmu1, mu1, mu2, var_as):

    alphaS = 0
    if var_as == 'nominal' : alphaS = cnst.asMZ
    elif var_as == 'up' : alphaS = cnst.asMZ+cnst.err_as
    elif var_as == 'down' : alphaS = cnst.asMZ-cnst.err_as
    else :
        print 'ERROR!'
        return 0

    crd = rundec.CRunDec()
    
    asmu1 = crd.AlphasExact(alphaS, cnst.MZ, mu1, cnst.nflav, cnst.nloops)
    asmu2 = crd.AlphasExact(alphaS, cnst.MZ, mu2, cnst.nflav, cnst.nloops)
    
    mtmu2 = crd.mMS2mMS(mtmu1, asmu1, asmu2, cnst.nflav, cnst.nloops)
    
    return mtmu2


################################

# "formInputFileName" provides the correct name of the .dat input file from where
# the calculated cross section is read out

################################


def formInputFileName ( renscale, facscale, topmass, pdfmember ):

    infileName='tt_tot_tot_ABMP16_'
    infileName+=str(pdfmember)+'_'
    if pdfmember<10 : infileName+='_'
    
    if renscale>=1000 : infileName+='_'
    if renscale>=100 : infileName+=str(renscale)+'_'
    else : 
        renstr = TString(str(renscale))
        renstr.ReplaceAll('.','_.')
        infileName+=str(renstr)+'_'

    if facscale>=1000 : infileName+='_'
    if facscale>=100 : infileName+=str(facscale)+'_'
    else : 
        facstr = TString(str(facscale))
        facstr.ReplaceAll('.','_.')
        infileName+=str(facstr)+'_'

    infileName+=str(topmass)+'_MSbar.dat'

    #because of some problem with MCFM outputs...
    if TString(infileName).Contains('.2') :
        infileName=str(TString(infileName).ReplaceAll('.2','.1'))
    if TString(infileName).Contains('.6') :
        infileName=str(TString(infileName).ReplaceAll('.6','.5'))

        
    return infileName


################################

# "readCalculatedXsec" reads the calculated cross section in a given mtt bin
# for any value of mu_r, mu_f, mt and PDF member
# calls the function "formInputFileName"

################################


def readCalculatedXsec (renscale, facscale, topmass, pdfmember, mttbin):

    fileName = formInputFileName ( renscale, facscale, topmass, pdfmember )
    if not os.path.isfile('out_hists/'+fileName):
        print 'WARNING: missing file', fileName
        return 0
    infile = open('out_hists/'+fileName,'r')
    lines = infile.read().splitlines()
    fillList = False
    bincenters = []
    xsecs = []
    for line in lines: 
        if 'Cross-section' in line:
            xsec_tot = line.split(':')
            xsec_tot = xsec_tot[1].split('+/-')
            xsec_tot = float(xsec_tot[0])/1000.
        if 'HIST =   2' in line: fillList = True
        if fillList and not 'm34' in line and not 'HIST' in line:
            bin_center, bin_xsec, bin_err = line.split()
            bin_center = float(bin_center)
            bin_xsec = float(bin_xsec)
            if bin_xsec == 0 and bin_center > 300 : fillList = False
            bincenters.append(bin_center)
            xsecs.append(bin_xsec)

    bin_width=bincenters[1]-bincenters[0]

    if mttbin < 3:
        xsec = xsecs[mttbin]/1000.*bin_width
    elif cnst.n_mttbins == 3:
        # xsec=xsec_tot-(xsecs[1]+xsecs[2])/1000.*bin_width
        xsec=0
        for i in range(mttbin,mttbin+11):
            xsec += xsecs[i]/1000.*bin_width
    elif mttbin==3: # n_mttbins = 4
        xsec = (xsecs[3]+xsecs[4])/1000.*bin_width
    else:
        xsec=0
        for i in range(mttbin+1,mttbin+16):
            xsec += xsecs[i]/1000.*bin_width
            
        
    return xsec 


################################

# "getMassAndError" provides the extracted mass and its error in a given mtt bin
# and for a choice of mu_r, mu_f, pdf member and extrapolation uncertainty

################################

def getMassAndError(mttbin, murscale, mufscale, pdfmember, extrapol, contrib): 

    graph = TGraph()
    i=-1
    chi2_v = []
    mass_v_sel = []
    for mass in mass_v:
        i+=1
        mur = mass
        muf = mass
        if (murscale == 'up')   : mur*=2
        if (mufscale == 'up')   : muf*=2
        if (murscale == 'down') : mur*=.5
        if (mufscale == 'down') : muf*=.5
        xsec_th = readCalculatedXsec(mur, muf, mass, pdfmember, mttbin) 
        if xsec_th == 0 : continue # torm
        xsec_exp = xsec_err = 0.
        if mttbin == 1 : 
            xsec_exp = xsec_1 
            xsec_err = ((xsec_1/100.*(var.err_xsec_1_up+var.err_xsec_1_down)/2.)**2 + var.err_xsec_toys_1**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + var.extr_1_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_1_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + var.extr_1_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_1_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + var.contribs_1[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - var.contribs_1[abs(contrib)-1]/100.
                
        elif mttbin == 2 : 
            xsec_exp = xsec_2 
            xsec_err = ((xsec_2/100.*(var.err_xsec_2_up+var.err_xsec_2_down)/2.)**2 + var.err_xsec_toys_2**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + var.extr_2_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_2_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + var.extr_2_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_2_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + var.contribs_2[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - var.contribs_2[abs(contrib)-1]/100.

        else : # mttbin = 3
            xsec_exp = xsec_3 
            xsec_err = ((xsec_3/100.*(var.err_xsec_3_up+var.err_xsec_3_down)/2.)**2 + var.err_xsec_toys_3**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + var.extr_3_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_3_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + var.extr_3_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_3_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + var.contribs_3[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - var.contribs_3[abs(contrib)-1]/100.


        # chi2 = abs(xsec_th-xsec_exp)/xsec_err
        chi2 = (xsec_th-xsec_exp)/xsec_err
        fact_A = 0
        if mttbin == 1 : fact_A = (var.err_xsec_1_up-var.err_xsec_1_down)/(var.err_xsec_1_up+var.err_xsec_1_down)
        elif mttbin == 2 : fact_A = (var.err_xsec_2_up-var.err_xsec_2_down)/(var.err_xsec_2_up+var.err_xsec_2_down)
        else : fact_A = (var.err_xsec_3_up-var.err_xsec_3_down)/(var.err_xsec_3_up+var.err_xsec_3_down)
        chi2 *= (1-2*fact_A*chi2+5*fact_A*fact_A*chi2*chi2)**.5
        # graph.SetPoint(i,mass,chi2)
        if mass != mass_v[0]:
            # if chi2<chi2_v[len(chi2_v)-1] : continue
            chi2_slope = (chi2-chi2_v[len(chi2_v)-1])/(mass-mass_v_sel[len(mass_v_sel)-1])
            slope_low = 0
            slope_high = 0
            if mttbin == 1 :
                slope_low = -6.0
                slope_high = -0.7
            elif mttbin == 2 :
                slope_low = -1.0
                slope_high = 0.5
            else:
                slope_low = -1.0
                slope_high = 0.05

            if chi2_slope < slope_low or chi2_slope > slope_high :
                if mass!= mass_v[1]: continue
                else:
                    if chi2 < chi2_v[len(chi2_v)-1] : continue

        if mttbin != 3:
            if mass < cnst.mass_fine_min : break
        mass_v_sel.append(mass)
        chi2_v.append(chi2)
            
    for i in range(0,len(chi2_v)):
        graph.SetPoint(i,mass_v_sel[i],chi2_v[i])

    funct = TF1('funct','pol4(0)',cnst.mass_min,cnst.mass_max)
    funct.SetParameter(0,1080)
    funct.SetParameter(1,-12)
    funct.SetParameter(2,0.033)
    # graph.Fit(funct,'q','',163,166)
    if mttbin==3: graph.Fit(funct,'q','',cnst.mass_min+0.1,cnst.mass_max-0.1)
    else: graph.Fit(funct,'q','',cnst.mass_fine_min+0.1,cnst.mass_max-0.1)
    
    if mttbin==3:
        line_up = TLine(cnst.mass_min,1.,cnst.mass_max,1.)
        line_down = TLine(cnst.mass_min,-1.,cnst.mass_max,-1.)
        line_central = TLine(cnst.mass_min,0.,cnst.mass_max,0.)
    else :
        line_up = TLine(cnst.mass_fine_min,1.,cnst.mass_fine_max,1.)
        line_down = TLine(cnst.mass_fine_min,-1.,cnst.mass_fine_max,-1.)
        line_central = TLine(cnst.mass_fine_min,0.,cnst.mass_fine_max,0.)


    line_up.SetLineColor(rt.kGreen+2)
    line_down.SetLineColor(rt.kGreen+2)
    line_central.SetLineColor(rt.kRed)

    outdir = 'plots_chi2'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c = TCanvas()
    c.SetGrid()
    graph.Draw('ap')
    graph.SetMarkerStyle(7)
    graph.SetTitle('signed #chi^{2} vs top mass in MCFM prediction; m_{t}(m_{t}) [GeV]; data-theory #chi^{2}')
    line_up.Draw('same')
    line_down.Draw('same')
    line_central.Draw('same')
    if not doingToys:
        c.Print(outdir+'/test_mtt'+str(mttbin)+'_mur_'+murscale+'_muf_'+mufscale+'_extrapol'+str(extrapol)+'_pdf'+str(pdfmember)+'_contrib'+str(contrib)+'.png')

    if mttbin==3:
        fitted_mass = funct.GetX(0,cnst.mass_min-30,cnst.mass_max+30)
        fitted_mass_up = funct.GetX(-1,fitted_mass-10,fitted_mass+10)
        fitted_mass_down = funct.GetX(1,fitted_mass-10,fitted_mass+10)
    elif mttbin==2:
        fitted_mass = funct.GetX(0,cnst.mass_fine_min,cnst.mass_fine_max)
        fitted_mass_up = funct.GetX(-1,fitted_mass-7,fitted_mass+7)
        fitted_mass_down = funct.GetX(1,fitted_mass-7,fitted_mass+7)
    else:        
        fitted_mass = funct.GetX(0,cnst.mass_fine_min,cnst.mass_fine_max)
        fitted_mass_up = funct.GetX(-1,fitted_mass-3,fitted_mass+3)
        fitted_mass_down = funct.GetX(1,fitted_mass-3,fitted_mass+3)
    #now evolve masses
    mu = 0
    if mttbin == 1 : mu = cnst.mu_1
    elif mttbin == 2 : mu = cnst.mu_2
    else : mu = cnst.mu_3 # mttbin = 3

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0 and contrib==0 and not doingToys:
        print
        print 'mt(mt) bin', mttbin,'=', round(fitted_mass,2), round(fitted_mass_up-fitted_mass,2), round(fitted_mass-fitted_mass_down,2)

    fitted_mtmt = fitted_mass
    fitted_mass = mtmt2mtmu(fitted_mass, mu)
    fitted_mass_up = mtmt2mtmu(fitted_mass_up, mu)
    fitted_mass_down = mtmt2mtmu(fitted_mass_down, mu)

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0 and contrib==0 and not doingToys:
        print 'mt(mu) bin', mttbin,'=', round(fitted_mass,2), round(fitted_mass_up-fitted_mass,2), round(fitted_mass-fitted_mass_down,2)
        
    fitted_mass_err = (fitted_mass_up - fitted_mass_down)/2

    return [fitted_mass, fitted_mass_err, fitted_mtmt]


################################

# "getRatios" calculates all the ratios and their experimental errors
# taking the correlations into account

################################

def getRatios(mass_1, mass_2, mass_3, err_1, err_2, err_3):

    ratio_1_2 = mass_1/mass_2
    err_ratio_1_2 = (err_1/mass_1)**2 + (err_2/mass_2)**2 - 2*var.corr_1_2*(err_1/mass_1)*(err_2/mass_2)
    err_ratio_1_2 = err_ratio_1_2**.5
    err_ratio_1_2*= ratio_1_2

    ratio_3_2 = mass_3/mass_2
    err_ratio_3_2 = (err_3/mass_3)**2 + (err_2/mass_2)**2 - 2*var.corr_3_2*(err_3/mass_3)*(err_2/mass_2)
    err_ratio_3_2 = err_ratio_3_2**.5
    err_ratio_3_2*= ratio_3_2

    return [ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2]
    

################################

# "getScaleUncertainties" calculates the scale uncertainties of the ratios

################################


def getScaleUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_scale_up_1_2 = central_ratio_1_2
    err_scale_up_3_2 = central_ratio_3_2
    err_scale_down_1_2 = central_ratio_1_2
    err_scale_down_3_2 = central_ratio_3_2

    for renscale in 'nominal', 'up','down':
        for facscale in 'nominal', 'up','down':

            if renscale == 'up' and facscale == 'down': continue
            if facscale == 'up' and renscale == 'down': continue
            if facscale == 'nominal' and renscale == 'nominal': continue
            
            mass_and_err_1 = getMassAndError(1, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_2 = getMassAndError(2, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_3 = getMassAndError(3, renscale, facscale, 0 , 0 , 0 )

            ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                        mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

            ratio_1_2 = ratios_and_errs[0]
            ratio_3_2 = ratios_and_errs[1]
    
            if ratio_1_2 > err_scale_up_1_2 : err_scale_up_1_2 = ratio_1_2
            if ratio_3_2 > err_scale_up_3_2 : err_scale_up_3_2 = ratio_3_2
            if ratio_1_2 < err_scale_down_1_2 : err_scale_down_1_2 = ratio_1_2
            if ratio_3_2 < err_scale_down_3_2 : err_scale_down_3_2 = ratio_3_2
    

    err_scale_up_1_2 -= central_ratio_1_2
    err_scale_up_3_2 -= central_ratio_3_2
    err_scale_down_1_2 = central_ratio_1_2 - err_scale_down_1_2
    err_scale_down_3_2 = central_ratio_3_2 - err_scale_down_3_2

    return [err_scale_up_1_2, err_scale_down_1_2, err_scale_up_3_2, err_scale_down_3_2]
    

################################

# "getPDFUncertainties" calculates the PDF uncertainties of the ratios

################################

def getPDFUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_pdf_up_1_2 = 0
    err_pdf_up_3_2 = 0
    err_pdf_down_1_2 = 0
    err_pdf_down_3_2 = 0
    
    for pdf in range(1,30):
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', pdf , 0 , 0 )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]

        err_pdf_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_pdf_3_2 = (ratio_3_2-central_ratio_3_2)**2

        if ratio_1_2 > central_ratio_1_2 : err_pdf_up_1_2 += err_pdf_1_2
        else : err_pdf_down_1_2 += err_pdf_1_2

        if ratio_3_2 > central_ratio_3_2 : err_pdf_up_3_2 += err_pdf_3_2
        else : err_pdf_down_3_2 += err_pdf_3_2
        
    return [err_pdf_up_1_2**.5, err_pdf_down_1_2**.5, err_pdf_up_3_2**.5, err_pdf_down_3_2**.5]
        

################################

# "getExtrapolationUncertainties" calculates the extrapolation uncertainties of the ratios

################################

def getExtrapolationUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_extr_up_1_2 = 0
    err_extr_up_3_2 = 0
    err_extr_down_1_2 = 0
    err_extr_down_3_2 = 0

    print '\n'
    print 'extrapol', 'ratio_1_2', 'ratio_3_2\n'
    for extr in range(-len(var.extr_1_up),len(var.extr_1_up)+1) :
        if extr == 0 : continue

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , extr , 0 )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]

        err_extr_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_extr_3_2 = (ratio_3_2-central_ratio_3_2)**2

        name = var.extr_name[abs(extr)-1]
        if extr>0 : name+='_up'
        else : name+='_down'
        print name, round(100*(ratio_1_2/central_ratio_1_2-1),2), round(100*(ratio_3_2/central_ratio_3_2-1),2),'%'

        
        if ratio_1_2 > central_ratio_1_2 : err_extr_up_1_2 += err_extr_1_2
        else : err_extr_down_1_2 += err_extr_1_2

        if ratio_3_2 > central_ratio_3_2 : err_extr_up_3_2 += err_extr_3_2
        else : err_extr_down_3_2 += err_extr_3_2


    return [err_extr_up_1_2**.5, err_extr_down_1_2**.5, err_extr_up_3_2**.5, err_extr_down_3_2**.5]


################################

# "makeTheoryPrediction" calculates the running of mt(mu) starting from m(mu_2)

################################


def makeTheoryPrediction(outfile, mass_2):

    scales = []
    r_up = []
    r_down = []
    out = open(outfile+'.txt','w')
        
    for scale in range(350,670+1):
        ratio = mtmu2mtmu(mass_2, cnst.mu_2, scale, 'nominal')/mass_2
        out.write(str(scale)+'\t'+str(ratio)+'\n')
        ratio_up = mtmu2mtmu(mass_2, cnst.mu_2, scale, 'up')/mass_2
        ratio_down = mtmu2mtmu(mass_2, cnst.mu_2, scale, 'down')/mass_2
        r_up.append(ratio_up)
        r_down.append(ratio_down)
        scales.append(scale)
        
    out.close()
        
    return [scales, r_up, r_down]



################################

# "makeAdditionalTheoryPrediction" calculates the running of mt(mu) starting from mt(mt)

################################


def makeAdditionalTheoryPrediction (mtmt, err_mtmt, mtmu):
    r = []
    ru = []
    rd = []
    scales = []
    r.append(mtmt/mtmu)
    ru.append((mtmt+err_mtmt)/mtmu)
    rd.append((mtmt-err_mtmt)/mtmu)
    scales.append(mtmt)
    for scale in range(int(mtmt)+1,670+1):
        ratio = mtmt2mtmu(mtmt,scale)/mtmu
        ratio_up = mtmt2mtmu(mtmt+err_mtmt,scale)/mtmu
        ratio_down = mtmt2mtmu(mtmt-err_mtmt,scale)/mtmu
        r.append(ratio)
        ru.append(ratio_up)
        rd.append(ratio_down)
        scales.append(scale)
        
    return [r, ru, rd, scales]



################################

# "makeRatioPlots" produces the ratio plots for the running of mt

################################


def makeRatioPlots (mass_2, ratio_12, ratio_32, err_12, err_32, mtmt_2):

    graph = TGraphErrors(3)
    
    graph.SetPoint(0,cnst.mu_1,ratio_12)
    graph.SetPointError(0,0,err_12)

    graph.SetPoint(1,cnst.mu_2,1)
    graph.SetPointError(1,0,0)

    graph.SetPoint(2,cnst.mu_3,ratio_32)
    graph.SetPointError(2,0,err_32)

    theoryFileName = 'theory_prediction'
    l = makeTheoryPrediction(theoryFileName,mass_2)
    th = TGraph(theoryFileName+'.txt')    
    th.SetLineColor(rt.kRed);
    th.SetLineWidth(2);

    os.remove(theoryFileName+'.txt')

    scales = l[0]
    r_up = l[1]
    r_down = l[2]
    
    th_band = TGraph (2*th.GetN())

    for i in range(0, th.GetN()):
        th_band.SetPoint(i,scales[i],max(r_up[i],r_down[i]))
        th_band.SetPoint(th.GetN()+i,scales[th.GetN()-i-1],min(r_up[th.GetN()-i-1],r_down[th.GetN()-i-1]));

    th_band.SetFillStyle(3001);
    th_band.SetFillColor(rt.kRed);
    
    leg = TLegend(.15,.2,.6,.32)
    leg.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg.AddEntry(th,'RunDec @ 2 loops (5 flav.)','l')

    latexLabel1 = TLatex();
    latexLabel1.SetTextSize(0.04);
    latexLabel1.SetTextFont(42)
    latexLabel1.SetNDC();
    
    graph.GetXaxis().SetTitle('centre-of-gravity of m_{t#bar{t}} [GeV]')
    graph.GetYaxis().SetTitle('m_{t}(m_{t#bar{t}}) / m_{t}('+str(int(cnst.mu_2))+' GeV)')
    graph.SetTitle('running of m_{t}(#mu) as a function of #mu=m_{t#bar{t}}')
    graph.SetMarkerStyle(8)
    
    g = graph.Clone()
    g.RemovePoint(1)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c = TCanvas()
    g.SetMarkerStyle(8)
    g.Draw('ap')
    th.Draw("L same");
    g1.Draw('p same')
    leg.Draw('same')
    th_band.Draw('f same')
    latexLabel1.DrawLatex(0.59, 0.78, "ABMP16_5_nlo PDF set");

    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/running.png')
    c.SaveAs(outdir+'/running.pdf')
    c.SaveAs(outdir+'/running.root')

    # mtmu = mtmt2mtmu (cnst.mtmt, cnst.mu_2)
    mtmu = mass_2

    gr_add = TGraphErrors()

    gr_add.SetPoint(0,cnst.mtmt,cnst.mtmt/mtmu)
    gr_add.SetPointError(0,0,cnst.mtmt_err/mtmu)

    l = makeAdditionalTheoryPrediction(cnst.mtmt, cnst.mtmt_err, mtmu)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]

    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1]);

    gr_band.SetFillStyle(3002);
    gr_band.SetFillColor(rt.kRed+1);

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(#mu) / m_{t}('+str(int(cnst.mu_2))+' GeV)')
    gr_band.SetTitle('running of m_{t}(#mu) as a function of #mu')

    gr_add.SetMarkerStyle(8)
    gr_add.SetMarkerColor(rt.kBlue)
    gr_add.SetLineColor(rt.kBlue)

    gr_band.GetYaxis().SetRangeUser(0.9,1.13)

    leg2 = leg.Clone()
    leg2.Clear()
    leg2.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg2.AddEntry(gr_add,'Hathor @NLO from incl. #sigma_{t#bar{t}} (same data)')
    leg2.AddEntry(gr_band,'RunDec @ 2 loops (5 flav.)','f')

    g = graph.Clone()
    g.RemovePoint(1)
    g.SetMarkerStyle(8)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c.Clear()
    gr_band.Draw('af')
    gr_add.Draw('p same')
    g.Draw('p same')
    g1.Draw('p same')
    leg2.Draw('same')

    c.SaveAs(outdir+'/test_incl.png')
    c.SaveAs(outdir+'/test_incl.pdf')
    c.SaveAs(outdir+'/test_incl.root')


    graph.Clear()

    graph.SetPoint(0,cnst.mu_1,ratio_12/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_1, 'nominal')/mass_2))
    graph.SetPointError(0,0,err_12/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_1, 'nominal')/mass_2))

    graph.SetPoint(1,cnst.mu_2,1)
    graph.SetPointError(1,0,0)

    graph.SetPoint(2,cnst.mu_3,ratio_32/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_3, 'nominal')/mass_2))
    graph.SetPointError(2,0,err_32/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_3, 'nominal')/mass_2))

    line = TLine(370,1,660,1)
    line.SetLineColor(rt.kRed)
    line.SetLineWidth(2)
    
    leg = TLegend(.15,.2,.7,.32)
    leg.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')

    graph.GetXaxis().SetTitle('centre-of-gravity of m_{t#bar{t}} [GeV]')
    graph.GetYaxis().SetTitle('m_{t}^{#mu}(m_{t}) / m_{t}^{#mu_{2}}(m_{t})')
    graph.SetTitle('extracted m_{t}(m_{t}) in bins of #mu=m_{t#bar{t}}')

    g = graph.Clone()
    g.RemovePoint(1)
    g.SetMarkerStyle(8)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c = TCanvas()
    g.Draw('ap')
    leg.Draw('same')
    line.Draw('same')
    g1.Draw('p same')
    latexLabel1.DrawLatex(0.59, 0.8, "ABMP16_5_nlo PDF set");

    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/mt_mt.png')
    c.SaveAs(outdir+'/mt_mt.pdf')
    c.SaveAs(outdir+'/mt_mt.root')

    
    gr_add.Clear()

    gr_add.SetPoint(0,cnst.mtmt,cnst.mtmt/mtmt_2)
    gr_add.SetPointError(0,0,cnst.mtmt_err/mtmt_2)

    l = makeAdditionalTheoryPrediction(cnst.mtmt, cnst.mtmt_err, mtmu)
    scales = l[3]

    gr_band.Clear()
    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],(cnst.mtmt+cnst.mtmt_err)/mtmt_2)
        gr_band.SetPoint(len(scales)+i,scales[len(scales)-i-1],(cnst.mtmt-cnst.mtmt_err)/mtmt_2);

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(m_{t}) / m_{t}^{#mu_{2}}(m_{t})')
    gr_band.SetTitle('m_{t}(m_{t}) measured as a function of #mu')

    gr_band.SetFillStyle(3002);
    gr_band.SetFillColor(rt.kRed+1);
    gr_band.GetYaxis().SetRangeUser(0.94,1.04)


    leg3 = TLegend(.15,.2,.7,.32)
    leg3.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg3.AddEntry(gr_add,'Hathor @NLO from incl. #sigma_{t#bar{t}} (same data)')
    leg3.AddEntry(gr_band,'experimental + extrapolation + PDF uncertainties','f')

    g = graph.Clone()
    g.RemovePoint(1)
    g.SetMarkerStyle(8)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c.Clear()
    gr_band.Draw('af')
    gr_add.Draw('p same')
    g.Draw('p same')
    g1.Draw('p same')
    leg3.Draw('same')

    c.SaveAs(outdir+'/test_mtmt.png')
    c.SaveAs(outdir+'/test_mtmt.pdf')
    c.SaveAs(outdir+'/test_mtmt.root')

    return



################################

# "makeChi2Test" calculates the significance of the obseved running
#  using a chi2 test

################################


def makeChi2Test(mass2, ratio12, ratio32, err12, err32):

    th_ratio12 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_1, 'nominal')/mass2
    th_ratio32 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_3, 'nominal')/mass2

    print 
    print 'theory:', round(th_ratio12,3), round(th_ratio32,3)
    
    chi2_no_run = ((ratio12-1)/err12)**2 + ((ratio32-1)/err32)**2
    chi2_run = ((ratio12-th_ratio12)/err12)**2 + ((ratio32-th_ratio32)/err32)**2

    print 
    print 'chi2: running hp', round(chi2_run,2)
    print 'chi2: no running hp', round(chi2_no_run,2)
    print
    print 'significance =', round((chi2_no_run-chi2_run)**.5,2)
    print
    
    
    return



################################

# "estimateSystContributions" estimates the contribution to the total uncertainty
#  of each individual systematics source

################################


def estimateSystContributions(central_ratio_1_2, central_ratio_3_2):
    
    print
    print 'contribs\n'

    for contrib in range(1,len(var.contribs_1)+1) :

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , contrib )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , contrib )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , contrib )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2_up = ratios_and_errs[0]
        ratio_3_2_up = ratios_and_errs[1]

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2_down = ratios_and_errs[0]
        ratio_3_2_down = ratios_and_errs[1]

        contr_1_2_up = (ratio_1_2_up/central_ratio_1_2-1)*100
        contr_1_2_down = (ratio_1_2_down/central_ratio_1_2-1)*100
        contr_3_2_up = (ratio_3_2_up/central_ratio_3_2-1)*100
        contr_3_2_down = (ratio_3_2_down/central_ratio_3_2-1)*100

        
        if (contr_1_2_up*contr_1_2_down) < 0:
            contr_1_2 = abs(contr_1_2_up-contr_1_2_down)/2.
        else : contr_1_2 = max(abs(contr_1_2_up),abs(contr_1_2_down))

        if (contr_3_2_up*contr_3_2_down) < 0:
            contr_3_2 = abs(contr_3_2_up-contr_3_2_down)/2.
        else : contr_3_2 = max(abs(contr_3_2_up),abs(contr_3_2_down))

        if contr_1_2 > 0.3 or contr_3_2 > 0.5 :
            print var.syst_names[contrib-1], round(contr_1_2,2), round(contr_3_2,2)


    print
    
    return



################################

# "throwToyCrossSections" throws toy cross sections taking the correlations into account
#  the results are used to estimate the correlation between the fitted masses

################################

def throwToyCrossSections(r):

    global xsec_1, xsec_2, xsec_3

    err_1 = xsec_1*(var.err_xsec_1_up+var.err_xsec_1_down)/2./100.
    err_2 = xsec_2*(var.err_xsec_2_up+var.err_xsec_2_down)/2./100.
    err_3 = xsec_3*(var.err_xsec_3_up+var.err_xsec_3_down)/2./100.

    V = [[err_1*err_1, var.corr_1_2*err_1*err_2, var.corr_1_3*err_1*err_3],
         [var.corr_1_2*err_1*err_2, err_2*err_2, var.corr_3_2*err_2*err_3],
         [var.corr_1_3*err_1*err_3, var.corr_3_2*err_2*err_3, err_3*err_3]]



    L = np.linalg.cholesky(V) # Cholesky decomposition 

    mu = [xsec_1, xsec_2, xsec_3]
    z = [r.Gaus(0,1),r.Gaus(0,1),r.Gaus(0,1)]

    y = L.dot(z)
    y += mu

    xsec_1 = y[0]
    xsec_2 = y[1]
    xsec_3 = y[2]
        
    return



################################

# "estimateMassCorrelations" estimates the correlation between the fitted masses
#  to run, set "ntoys" to some positive numbers (at least 10k)
#  N.B. it is absolutely safe to use the correlaton between the cross section (i.e. ntoys=0)

################################


def estimateMassCorrelations():

    global doingToys
    global xsec_1, xsec_2, xsec_3
    
    orig_xsec_1 = xsec_1
    orig_xsec_2 = xsec_2
    orig_xsec_3 = xsec_3
    
    doingToys = True
    
    r=TRandom3()

    m12 = TH2D('m12','m12',100,150,160,100,130,170)
    m32 = TH2D('m32','m32',100,120,180,100,130,170)

    xs12 = TH2D('xs12','xs12',100,xsec_1*(1-5*var.err_xsec_1_down/100.),xsec_1*(1+5*var.err_xsec_1_up/100.),
                100,xsec_2*(1-5*var.err_xsec_2_down/100.),xsec_2*(1-5*var.err_xsec_2_down/100.))
    xs13 = TH2D('xs13','xs13',100,xsec_1*(1-5*var.err_xsec_1_down/100.),xsec_1*(1+5*var.err_xsec_1_up/100.),
                100,xsec_3*(1-5*var.err_xsec_3_down/100.),xsec_3*(1-5*var.err_xsec_3_down/100.))
    xs23 = TH2D('xs23','xs23',100,xsec_2*(1-5*var.err_xsec_2_down/100.),xsec_2*(1+5*var.err_xsec_2_up/100.),
                100,xsec_3*(1-5*var.err_xsec_3_down/100.),xsec_3*(1-5*var.err_xsec_3_down/100.))


    
    r_corr = TH2D('r_corr','r_corr',100,0.9,1.15,100,0.75,1.15)
    
    
    for i in range(1,ntoys+1):
        if i%1000==0 or i==1: print 'toy n.', i
        throwToyCrossSections(r)

        xs12.Fill(xsec_1,xsec_2)
        xs13.Fill(xsec_1,xsec_3)
        xs23.Fill(xsec_2,xsec_3)
        
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , 0)
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , 0)
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , 0)
        m12.Fill(mass_and_err_1[0],mass_and_err_2[0])
        m32.Fill(mass_and_err_3[0],mass_and_err_2[0])

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        r_corr.Fill(ratios_and_errs[0],ratios_and_errs[1])
        
        #extremely important
        xsec_1 = orig_xsec_1
        xsec_2 = orig_xsec_2
        xsec_3 = orig_xsec_3
    
    doingToys = False

    if replace_corr:
        var.corr_1_2 = m12.GetCorrelationFactor()
        var.corr_3_2 = m32.GetCorrelationFactor()
    
    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c = TCanvas()
    m12.SetTitle('correlation = '+str(round(m12.GetCorrelationFactor(),2)))
    m12.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_1))+' GeV) [GeV]')
    m12.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_2))+' GeV) [GeV]')
    m12.DrawNormalized('colz')
    c.SaveAs(outdir+'/mass_corr_1_2.png')
    c.SaveAs(outdir+'/mass_corr_1_2.pdf')
    c.SaveAs(outdir+'/mass_corr_1_2.root')
    c.Clear()
    m32.SetTitle('correlation = '+str(round(m32.GetCorrelationFactor(),2)))
    m32.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_3))+' GeV) [GeV]')
    m32.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_2))+' GeV) [GeV]')
    m32.DrawNormalized('colz')
    c.SaveAs(outdir+'/mass_corr_3_2.png')
    c.SaveAs(outdir+'/mass_corr_3_2.pdf')
    c.SaveAs(outdir+'/mass_corr_3_2.root')
    c.Clear()
    r_corr.SetTitle('correlation = '+str(round(r_corr.GetCorrelationFactor(),2)))
    r_corr.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_1))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_3))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr.DrawNormalized('colz')
    c.SaveAs(outdir+'/ratio_corr.png')
    c.SaveAs(outdir+'/ratio_corr.pdf')
    c.SaveAs(outdir+'/ratio_corr.root')

    xs12.SetTitle('correlation = '+str(round(xs12.GetCorrelationFactor(),2)))
    xs12.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{1}) [pb]')
    xs12.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{2}) [pb]')
    xs12.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_12.png')
    c.SaveAs(outdir+'/xsec_corr_12.pdf')
    c.SaveAs(outdir+'/xsec_corr_12.root')
    xs13.SetTitle('correlation = '+str(round(xs13.GetCorrelationFactor(),2)))
    xs13.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{1}) [pb]')
    xs13.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{3}) [pb]')
    xs13.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_13.png')
    c.SaveAs(outdir+'/xsec_corr_13.pdf')
    c.SaveAs(outdir+'/xsec_corr_13.root')
    xs23.SetTitle('correlation = '+str(round(xs23.GetCorrelationFactor(),2)))
    xs23.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{2}) [pb]')
    xs23.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{3}) [pb]')
    xs23.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_23.png')
    c.SaveAs(outdir+'/xsec_corr_23.pdf')
    c.SaveAs(outdir+'/xsec_corr_23.root')

    
    print
    if replace_corr:
        print 'new corr_1_2 =', round(var.corr_1_2,3)
        print 'new corr_3_2 =', round(var.corr_3_2,3)
    else:    
        print 'estimated corr_1_2 =', round(m12.GetCorrelationFactor(),3)
        print 'estimated corr_3_2 =', round(m32.GetCorrelationFactor(),3)
    print 'estimated corr ratios  =', round(r_corr.GetCorrelationFactor(),3)
    print
    
    return

################################

# main program

################################


def execute():

    global doingToys
    doingToys = False

    setMasses()
    if ntoys>0 : estimateMassCorrelations()

    mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , 0)
    mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , 0)
    mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , 0)

    ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

    ratio_1_2 = ratios_and_errs[0]
    ratio_3_2 = ratios_and_errs[1]
    err_ratio_1_2 = ratios_and_errs[2]
    err_ratio_3_2 = ratios_and_errs[3]


    pdf_errors = getPDFUncertainties(ratio_1_2, ratio_3_2)

    err_pdf_1_2_up = pdf_errors[0]
    err_pdf_1_2_down = pdf_errors[1]
    err_pdf_3_2_up = pdf_errors[2]
    err_pdf_3_2_down = pdf_errors[3]


    extr_errors = getExtrapolationUncertainties(ratio_1_2, ratio_3_2)

    err_extr_1_2_up = extr_errors[0]
    err_extr_1_2_down = extr_errors[1]
    err_extr_3_2_up = extr_errors[2]
    err_extr_3_2_down = extr_errors[3]

    err_1_2_up = (err_ratio_1_2**2 + err_pdf_1_2_up**2 + err_extr_1_2_up **2)**.5
    err_1_2_down = (err_ratio_1_2**2 + err_pdf_1_2_down**2 + err_extr_1_2_down **2)**.5

    err_3_2_up = (err_ratio_3_2**2 + err_pdf_3_2_up**2 + err_extr_3_2_up **2)**.5
    err_3_2_down = (err_ratio_3_2**2 + err_pdf_3_2_down**2 + err_extr_3_2_down **2)**.5

    print '\n'
    print 'uncertainties ratio_1_2:\n'
    print 'experimental =', round(err_ratio_1_2/ratio_1_2*100.,2), '%'
    print 'PDFs up/down =', round(err_pdf_1_2_up/ratio_1_2*100.,2), round(err_pdf_1_2_down/ratio_1_2*100.,2), '%'
    print 'extr up/down =', round(err_extr_1_2_up/ratio_1_2*100.,2), round(err_extr_1_2_down/ratio_1_2*100.,2), '%'
    print 'total =', round(max(err_1_2_up,err_1_2_down)/ratio_1_2*100.,2), '%'
    print '\n'
    print 'uncertainties ratio_3_2:\n'
    print 'experimental =', round(err_ratio_3_2/ratio_3_2*100.,2)
    print 'PDFs up/down =', round(err_pdf_3_2_up/ratio_3_2*100.,2), round(err_pdf_3_2_down/ratio_3_2*100.,2), '%'
    print 'extr up/down =', round(err_extr_3_2_up/ratio_3_2*100.,2), round(err_extr_3_2_down/ratio_3_2*100.,2), '%'
    print 'total =', round(max(err_3_2_up,err_3_2_down)/ratio_3_2*100.,2), '%'
    print '\n'

    print 'results:\n'
    print 'ratio_1_2 =', round(ratio_1_2,3), '+' , round(err_1_2_up,3), '-' , round(err_1_2_down,3)
    print 'ratio_3_2 =', round(ratio_3_2,3), '+' , round(err_3_2_up,3), '-' , round(err_3_2_down,3) 
    print

    makeRatioPlots (mass_and_err_2[0], ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2, mass_and_err_2[2])
    makeChi2Test (mass_and_err_2[0], ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2)
    if estimate_contribs : estimateSystContributions (ratio_1_2,ratio_3_2)
    
    return


################################

# execute main program

################################


execute()









