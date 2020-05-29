
import os
import rundec
import numpy as np

import ROOT as rt
import variables as var
import constants as cnst

from ROOT import TString, TH2D, TRandom3, TF1, TGraph, TLine, TCanvas, TGraphErrors, TGraphAsymmErrors, TLegend, TLatex
from variables import xsec_1, xsec_2, xsec_3, xsec_4


# options

estimate_contribs = False
estimate_significance = False

do_scale_variations = False # not yet well tested
facscale_only = True # only factorization scale

ntoys = 0
replace_corr = False  #recommended: False

preliminary = False

newscales = True

# end options



################################

# "setMasses" initializes which masses should be considered

################################

def setMasses():

    global mass_v
    mass_v = []
    if not newscales: i_mass = cnst.mass_max
    else: i_mass = cnst.mass_max_scales
    if not newscales: m_min = cnst.mass_min
    else: m_min = cnst.mass_min_scales
    while i_mass >= m_min :
        mass_v.append(i_mass)
        if i_mass <= cnst.mass_fine_max and i_mass > cnst.mass_fine_min+.1 and not newscales: i_mass-= 0.2
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

    topmassevolved=topmass
    if renscale != topmass:
        topmassevolved = float(int(mtmt2mtmu(topmass, renscale)*10))/10.

    origrenscale = renscale
    origfacscale = facscale
    
    if renscale != topmass or facscale != topmass:
        if TString(str(topmass)).EndsWith('.6') or TString(str(topmass)).EndsWith('.2'):
            renscale = round(renscale-.1,1)
            facscale = round(facscale-.1,1)
        elif TString(str(topmass)).EndsWith('.5'):
            if renscale < topmass : renscale = round(renscale-.1,1)
            if facscale < topmass : facscale = round(facscale-.1,1)
        
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

    tempname = infileName
    infileName+=str(topmassevolved)+'_MSbar.dat'

    #because of some problem with MCFM outputs...
    if renscale == facscale == topmass:
        if TString(infileName).Contains('.2') :
            infileName=str(TString(infileName).ReplaceAll('.2','.1'))
        if TString(infileName).Contains('.6') :
            infileName=str(TString(infileName).ReplaceAll('.6','.5'))
    elif origrenscale == topmass:
        if TString(infileName).Contains('.2_MS') :
            infileName=str(TString(infileName).ReplaceAll('.2_MS','.1_MS'))
        if TString(infileName).Contains('.6_MS') :
            infileName=str(TString(infileName).ReplaceAll('.6_MS','.5_MS'))
    else:
        if not os.path.isfile('out_scales/'+infileName):
            infileName=tempname+str(topmassevolved+.1)+'_MSbar.dat'
        if not os.path.isfile('out_scales/'+infileName):
            infileName=tempname+str(topmassevolved-.1)+'_MSbar.dat'
    return infileName


################################

# "newInputFileName" provides the correct name of the .dat input file from where
# the calculated cross section is read out

################################


def newInputFileName ( renscale, facscale, topmass, pdfmember ):

    infileName='tt_tot_tot_ABMP16_'
    infileName+=str(pdfmember)+'_'
    if pdfmember<10 : infileName+='_'

    topmassevolved=topmass
    if renscale != topmass:
        topmassevolved = float(int(mtmt2mtmu(topmass, renscale)*10))/10.

    renscale = round(renscale,0)
    facscale = round(facscale,0)
    
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

    tempname = infileName
    infileName+=str(topmassevolved)+'_MSbar.dat'

    if not os.path.isfile('scales_Grazzini/'+infileName):
        infileName=tempname+str(topmassevolved+.1)+'_MSbar.dat'
    if not os.path.isfile('scales_Grazzini/'+infileName):
        infileName=tempname+str(topmassevolved-.1)+'_MSbar.dat'

    return infileName

################################

# "readCalculatedXsec" reads the calculated cross section in a given mtt bin
# for any value of mu_r, mu_f, mt and PDF member
# calls the function "formInputFileName"

################################


def readCalculatedXsec (renscale, facscale, topmass, pdfmember, mttbin):

    if not newscales: fileName = formInputFileName ( renscale, facscale, topmass, pdfmember )
    else: fileName = newInputFileName ( renscale, facscale, topmass, pdfmember )
    if not newscales:
        if renscale == topmass and facscale == topmass: indir = 'out_hists/'
        else : indir = 'out_scales/'
    else: indir = 'scales_Grazzini/'
    if not os.path.isfile(indir+fileName):
        print 'WARNING: missing file', fileName
        return 0
    infile = open(indir+fileName,'r')
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
        if mttbin == 1: xsec += xsecs[0]/1000.*bin_width
    elif cnst.n_mttbins == 3:
        # xsec=xsec_tot-(xsecs[1]+xsecs[2])/1000.*bin_width
        xsec=0
        for i in range(mttbin,mttbin+11):
            xsec += xsecs[i]/1000.*bin_width
    elif mttbin==3: # n_mttbins = 4
        xsec = (xsecs[3]+xsecs[4])/1000.*bin_width
    else:
        xsec=xsec_tot
        for i in range(0,5):
            xsec -= xsecs[i]/1000.*bin_width
            
        
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
        if not newscales:
            mur = mass
            muf = mass
        else :
            if mttbin == 1:
                mur = cnst.mu_1
                muf = cnst.mu_1
            if mttbin == 2:
                mur = cnst.mu_2
                muf = cnst.mu_2
            if mttbin == 3:
                mur = cnst.mu_3
                muf = cnst.mu_3
            if mttbin == 4:
                mur = cnst.mu_4
                muf = cnst.mu_4
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

        elif mttbin == 3 : 
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

        elif mttbin == 4 : 
            xsec_exp = xsec_4 
            xsec_err = ((xsec_4/100.*(var.err_xsec_4_up+var.err_xsec_4_down)/2.)**2 + var.err_xsec_toys_4**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + var.extr_4_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_4_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + var.extr_4_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + var.extr_4_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + var.contribs_4[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - var.contribs_4[abs(contrib)-1]/100.


        # chi2 = abs(xsec_th-xsec_exp)/xsec_err
        chi2 = (xsec_th-xsec_exp)/xsec_err
        fact_A = 0
        if mttbin == 1 : fact_A = (var.err_xsec_1_up-var.err_xsec_1_down)/(var.err_xsec_1_up+var.err_xsec_1_down)
        elif mttbin == 2 : fact_A = (var.err_xsec_2_up-var.err_xsec_2_down)/(var.err_xsec_2_up+var.err_xsec_2_down)
        elif mttbin == 3 : fact_A = (var.err_xsec_3_up-var.err_xsec_3_down)/(var.err_xsec_3_up+var.err_xsec_3_down)
        elif mttbin == 4 : fact_A = (var.err_xsec_4_up-var.err_xsec_4_down)/(var.err_xsec_4_up+var.err_xsec_4_down)
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
            elif mttbin == 3 :
                slope_low = -1.0
                slope_high = 0.05
            elif mttbin == 4 :
                slope_low = -0.4
                slope_high = 0.2
                
            if chi2_slope < slope_low or chi2_slope > slope_high :
                if mass!= mass_v[1]: continue
                else:
                    if chi2 < chi2_v[len(chi2_v)-1] : continue

        if mttbin < 3 :
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
    if not newscales:
        if mttbin>=3: graph.Fit(funct,'q','',cnst.mass_min+0.1,cnst.mass_max-0.1)
        else: graph.Fit(funct,'q','',cnst.mass_fine_min+0.1,cnst.mass_max-0.1)
    else: graph.Fit(funct,'q','',cnst.mass_min_scales+0.5,cnst.mass_max_scales-0.5)
        
    if mttbin >= 3:
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
    elif mttbin==4:
        fitted_mass = funct.GetX(0,cnst.mass_min-10,cnst.mass_max+10)
        fitted_mass_up = funct.GetX(-1,fitted_mass-15,fitted_mass+15)
        fitted_mass_down = funct.GetX(1,fitted_mass-15,fitted_mass+15)
    else:        
        fitted_mass = funct.GetX(0,cnst.mass_fine_min,cnst.mass_fine_max)
        fitted_mass_up = funct.GetX(-1,fitted_mass-3,fitted_mass+3)
        fitted_mass_down = funct.GetX(1,fitted_mass-3,fitted_mass+3)
    #now evolve masses
    mu = 0
    if mttbin == 1 : mu = cnst.mu_1
    elif mttbin == 2 : mu = cnst.mu_2
    elif mttbin == 3 : mu = cnst.mu_3
    elif mttbin == 4 : mu = cnst.mu_4

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

def getRatios(mass_1, mass_2, mass_3, mass_4, err_1, err_2, err_3, err_4):

    ratio_1_2 = mass_1/mass_2
    err_ratio_1_2 = (err_1/mass_1)**2 + (err_2/mass_2)**2 - 2*var.corr_1_2*(err_1/mass_1)*(err_2/mass_2)
    err_ratio_1_2 = err_ratio_1_2**.5
    err_ratio_1_2*= ratio_1_2

    ratio_3_2 = mass_3/mass_2
    err_ratio_3_2 = (err_3/mass_3)**2 + (err_2/mass_2)**2 - 2*var.corr_3_2*(err_3/mass_3)*(err_2/mass_2)
    err_ratio_3_2 = err_ratio_3_2**.5
    err_ratio_3_2*= ratio_3_2

    ratio_4_2 = mass_4/mass_2
    err_ratio_4_2 = (err_4/mass_4)**2 + (err_2/mass_2)**2 - 2*var.corr_4_2*(err_4/mass_4)*(err_2/mass_2)
    err_ratio_4_2 = err_ratio_4_2**.5
    err_ratio_4_2*= ratio_4_2

    return [ratio_1_2, ratio_3_2, ratio_4_2, err_ratio_1_2, err_ratio_3_2, err_ratio_4_2]
    

################################

# "getScaleUncertainties" calculates the scale uncertainties of the ratios

################################


def getScaleUncertainties (central_ratio_1_2, central_ratio_3_2, central_ratio_4_2):

    err_scale_up_1_2 = central_ratio_1_2
    err_scale_up_3_2 = central_ratio_3_2
    err_scale_up_4_2 = central_ratio_4_2
    err_scale_down_1_2 = central_ratio_1_2
    err_scale_down_3_2 = central_ratio_3_2
    err_scale_down_4_2 = central_ratio_4_2

    for renscale in 'nominal', 'up','down':
        for facscale in 'nominal', 'up','down':

            if renscale == 'up' and facscale == 'down': continue
            if facscale == 'up' and renscale == 'down': continue
            if facscale == 'nominal' and renscale == 'nominal': continue

            if facscale_only:
                if renscale != 'nominal': continue # only factorization scale
            
            mass_and_err_1 = getMassAndError(1, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_2 = getMassAndError(2, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_3 = getMassAndError(3, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_4 = getMassAndError(4, renscale, facscale, 0 , 0 , 0 )

            ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                        mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

            ratio_1_2 = ratios_and_errs[0]
            ratio_3_2 = ratios_and_errs[1]
            ratio_4_2 = ratios_and_errs[2]
    
            if ratio_1_2 > err_scale_up_1_2 : err_scale_up_1_2 = ratio_1_2
            if ratio_3_2 > err_scale_up_3_2 : err_scale_up_3_2 = ratio_3_2
            if ratio_4_2 > err_scale_up_4_2 : err_scale_up_4_2 = ratio_4_2
            if ratio_1_2 < err_scale_down_1_2 : err_scale_down_1_2 = ratio_1_2
            if ratio_3_2 < err_scale_down_3_2 : err_scale_down_3_2 = ratio_3_2
            if ratio_4_2 < err_scale_down_4_2 : err_scale_down_4_2 = ratio_4_2
    

    err_scale_up_1_2 -= central_ratio_1_2
    err_scale_up_3_2 -= central_ratio_3_2
    err_scale_up_4_2 -= central_ratio_4_2
    err_scale_down_1_2 = central_ratio_1_2 - err_scale_down_1_2
    err_scale_down_3_2 = central_ratio_3_2 - err_scale_down_3_2
    err_scale_down_4_2 = central_ratio_4_2 - err_scale_down_4_2

    return [err_scale_up_1_2, err_scale_down_1_2, err_scale_up_3_2, err_scale_down_3_2, err_scale_up_4_2, err_scale_down_4_2]


################################

# "getPDFUncertainties" calculates the PDF uncertainties of the ratios

################################

def getPDFUncertainties (central_ratio_1_2, central_ratio_3_2, central_ratio_4_2):

    err_pdf_up_1_2 = 0
    err_pdf_up_3_2 = 0
    err_pdf_up_4_2 = 0
    err_pdf_down_1_2 = 0
    err_pdf_down_3_2 = 0
    err_pdf_down_4_2 = 0
    
    for pdf in range(1,30):
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', pdf , 0 , 0 )
        
        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]
        ratio_4_2 = ratios_and_errs[2]

        err_pdf_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_pdf_3_2 = (ratio_3_2-central_ratio_3_2)**2
        err_pdf_4_2 = (ratio_4_2-central_ratio_4_2)**2

        if ratio_1_2 > central_ratio_1_2 : err_pdf_up_1_2 += err_pdf_1_2
        else : err_pdf_down_1_2 += err_pdf_1_2

        if ratio_3_2 > central_ratio_3_2 : err_pdf_up_3_2 += err_pdf_3_2
        else : err_pdf_down_3_2 += err_pdf_3_2

        if ratio_4_2 > central_ratio_4_2 : err_pdf_up_4_2 += err_pdf_4_2
        else : err_pdf_down_4_2 += err_pdf_4_2

        
    return [err_pdf_up_1_2**.5, err_pdf_down_1_2**.5, err_pdf_up_3_2**.5, err_pdf_down_3_2**.5, err_pdf_up_4_2**.5, err_pdf_down_4_2**.5]
        

################################

# "getExtrapolationUncertainties" calculates the extrapolation uncertainties of the ratios

################################

def getExtrapolationUncertainties (central_ratio_1_2, central_ratio_3_2, central_ratio_4_2):

    err_extr_up_1_2 = 0
    err_extr_up_3_2 = 0
    err_extr_up_4_2 = 0
    err_extr_down_1_2 = 0
    err_extr_down_3_2 = 0
    err_extr_down_4_2 = 0

    print '\n'
    print 'extrapol', 'ratio_1_2', 'ratio_3_2', 'ratio_4_2\n'
    for extr in range(-len(var.extr_1_up),len(var.extr_1_up)+1) :
        if extr == 0 : continue

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , extr , 0 )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]
        ratio_4_2 = ratios_and_errs[2]

        err_extr_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_extr_3_2 = (ratio_3_2-central_ratio_3_2)**2
        err_extr_4_2 = (ratio_4_2-central_ratio_4_2)**2

        name = var.extr_name[abs(extr)-1]
        if extr>0 : name+='_up'
        else : name+='_down'
        print name, round(100*(ratio_1_2/central_ratio_1_2-1),2), round(100*(ratio_3_2/central_ratio_3_2-1),2), round(100*(ratio_4_2/central_ratio_4_2-1),2), '%'

        
        if ratio_1_2 > central_ratio_1_2 : err_extr_up_1_2 += err_extr_1_2
        else : err_extr_down_1_2 += err_extr_1_2

        if ratio_3_2 > central_ratio_3_2 : err_extr_up_3_2 += err_extr_3_2
        else : err_extr_down_3_2 += err_extr_3_2

        if ratio_4_2 > central_ratio_4_2 : err_extr_up_4_2 += err_extr_4_2
        else : err_extr_down_4_2 += err_extr_4_2


    return [err_extr_up_1_2**.5, err_extr_down_1_2**.5, err_extr_up_3_2**.5, err_extr_down_3_2**.5, err_extr_up_4_2**.5, err_extr_down_4_2**.5]


################################

# "makeTheoryPrediction" calculates the running of mt(mu) starting from m(mu_2)

################################


def makeTheoryPrediction(outfile, mass_2):

    scales = []
    r_up = []
    r_down = []
    out = open(outfile+'.txt','w')
        
    for scale in range(350/2,(1050+1)/2):
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


def makeAdditionalTheoryPrediction (mtmt, err_mtmt_up, err_mtmt_down, mtmu, doratio):
    r = []
    ru = []
    rd = []
    scales = []
    r.append(mtmt/mtmu)
    ru.append((mtmt+err_mtmt_up)/mtmu)
    rd.append((mtmt-err_mtmt_down)/mtmu)
    scales.append(mtmt)
    for scale in range(int(mtmt)+1,1050+1):
        ratio = mtmt2mtmu(mtmt,scale)/mtmu
        ratio_up = mtmt2mtmu(mtmt+err_mtmt_up,scale)/mtmu
        ratio_down = mtmt2mtmu(mtmt-err_mtmt_down,scale)/mtmu
        r.append(ratio)
        ru.append(ratio_up)
        rd.append(ratio_down)
        scales.append(scale)

    if not doratio:
        for i in range(0,len(r)):
            r[i]  *= mtmu
            ru[i] *= mtmu
            rd[i] *= mtmu
        
    return [r, ru, rd, scales]



################################

# "makeRatioPlots" produces the ratio plots for the running of mt

################################


def makeRatioPlots (mass_2, ratio_12, ratio_32, ratio_42, err_12_up, err_12_down, err_32_up, err_32_down, err_42_up, err_42_down, mtmt_2):

    graph = TGraphAsymmErrors(3)
    
    graph.SetPoint(0,cnst.mu_1,ratio_12)
    graph.SetPointError(0,0,0,err_12_down,err_12_up)

    graph.SetPoint(1,cnst.mu_2,1)
    graph.SetPointError(1,0,0,0,0)

    graph.SetPoint(2,cnst.mu_3,ratio_32)
    graph.SetPointError(2,0,0,err_32_down,err_32_up)

    graph.SetPoint(3,cnst.mu_4,ratio_42)
    graph.SetPointError(3,0,0,err_42_down,err_42_up)

    theoryFileName = 'theory_prediction'
    l = makeTheoryPrediction(theoryFileName,mass_2)
    th = TGraph(theoryFileName+'.txt')    
    th.SetLineColor(rt.kRed)
    th.SetLineWidth(2)

    os.remove(theoryFileName+'.txt')

    scales = l[0]
    r_up = l[1]
    r_down = l[2]
    
    th_band = TGraph (2*th.GetN())

    for i in range(0, th.GetN()):
        th_band.SetPoint(i,scales[i],max(r_up[i],r_down[i]))
        th_band.SetPoint(th.GetN()+i,scales[th.GetN()-i-1],min(r_up[th.GetN()-i-1],r_down[th.GetN()-i-1]))

    th_band.SetFillStyle(3001)
    th_band.SetFillColor(rt.kRed)
    

    latexLabel1 = TLatex()
    latexLabel1.SetTextSize(0.06)
    latexLabel1.SetNDC()
    
    latexLabel2 = TLatex()
    latexLabel2.SetTextSize(0.04)
    latexLabel2.SetTextFont(42)
    latexLabel2.SetNDC()
    
    graph.GetXaxis().SetTitle('#mu [GeV]')
    graph.GetXaxis().SetTitleSize(0.06)
    graph.GetXaxis().SetTitleOffset(.8)
    graph.GetYaxis().SetTitleSize(0.05)
    graph.GetYaxis().SetTitleOffset(1.2)
    graph.GetYaxis().SetLabelSize(0.05)
    graph.GetXaxis().SetLabelSize(0.05)
    graph.GetYaxis().SetTitle('m_{t}(#mu) / m_{t}(#mu_{ref})')
    graph.SetTitle('')
    graph.SetMarkerStyle(8)
    
    g = graph.Clone()
    g.RemovePoint(1)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    # g1.SetMarkerStyle(8)
    g1.SetMarkerStyle(4)
    # g1.SetMarkerSize(1.5)

    leg = TLegend(.17,.2,.78,.36)
    leg.SetBorderSize(0)
    leg.AddEntry(graph,'NLO extraction from differential #sigma_{t#bar{t}}','pe')
    leg.AddEntry(g1,'Reference scale (#mu = #mu_{ref})','p')
    leg.AddEntry(th,'one-loop RGE, n_{f} = 5, #alpha_{s}(m_{Z}) = 0.1191','l')
    
    c = TCanvas()
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.08)
    g.SetMarkerStyle(8)
    g.Draw('ap')
    th.Draw("L same")
    g1.Draw('p same')
    leg.Draw('same')
    th_band.Draw('f same')
    g.Draw('psame')
    latexLabel1.DrawLatex(0.16, 0.94, "CMS")
    latexLabel2.DrawLatex(0.76, 0.94, "35.9 fb^{-1} (13 TeV)")
    latexLabel2.DrawLatex(0.63, 0.81, "ABMP16_5_nlo PDF set")
    latexLabel2.DrawLatex(0.63, 0.76, "#mu_{ref} = "+str(int(cnst.mu_2))+" GeV")
    latexLabel2.DrawLatex(0.63, 0.71, "#mu_{0} = #mu_{ref}")
    if preliminary: latexLabel2.DrawLatex(0.205, 0.92 , "#it{Preliminary}")
    
    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/running.png')
    c.SaveAs(outdir+'/running.pdf')
    c.SaveAs(outdir+'/running.root')
    c.SaveAs(outdir+'/running.C')

    # mtmu = mtmt2mtmu (cnst.mtmt, cnst.mu_2)
    mtmu = mass_2

    gr_add = TGraphErrors()

    gr_add.SetPoint(0,cnst.mtmt,cnst.mtmt/mtmu)
    gr_add.SetPointError(0,0,cnst.mtmt_err/mtmu)

    l = makeAdditionalTheoryPrediction(cnst.mtmt, cnst.mtmt_err, cnst.mtmt_err, mtmu, True)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]

    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1])

    gr_band.SetFillStyle(3002)
    gr_band.SetFillColor(rt.kRed+1)

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(#mu) / m_{t}(#mu_{ref})')
    gr_band.SetTitle('')


    gr_band.GetXaxis().SetTitleSize(0.06)
    gr_band.GetXaxis().SetTitleOffset(.8)
    gr_band.GetYaxis().SetTitleSize(0.05)
    gr_band.GetYaxis().SetTitleOffset(1.1)
    gr_band.GetYaxis().SetLabelSize(0.05)
    gr_band.GetXaxis().SetLabelSize(0.05)

    
    gr_add.SetMarkerStyle(21)
    gr_add.SetMarkerColor(rt.kBlue)
    gr_add.SetLineColor(rt.kBlue)

    gr_band.GetYaxis().SetRangeUser(0.9,1.13)

    g = graph.Clone()
    g.RemovePoint(1)
    g1 = TGraph()
    g1.SetPoint(0,cnst.mu_2,1)
    g1.SetMarkerStyle(4)


    leg2 = TLegend(.14,.18,.86,.39)
    # leg2 = TLegend(.13,.15,.75,.32)
    leg2.SetBorderSize(0)
    leg2.AddEntry(graph,'NLO extraction from differential #sigma_{t#bar{t}}','pe')
    leg2.AddEntry(g1,'Reference value (#mu = #mu_{ref})','p')
    leg2.AddEntry(gr_add,'NLO extraction from inclusive #sigma_{t#bar{t}} (same data)','pe')
    leg2.AddEntry(gr_band,'Evolved uncert. One-loop RGE, n_{f} = 5, #alpha_{s}(m_{Z}) = 0.1191','f')


    
    c.Clear()
    gr_band.SetMinimum(0.79)
    gr_band.Draw('af')
    gr_add.Draw('p same')
    g.Draw('p same')
    g1.Draw('p same')
    leg2.Draw('same')
    latexLabel1.DrawLatex(0.16, 0.94, "CMS")
    latexLabel2.DrawLatex(0.76, 0.94, "35.9 fb^{-1} (13 TeV)")
    latexLabel2.DrawLatex(0.63, 0.81, "ABMP16_5_nlo PDF set")
    latexLabel2.DrawLatex(0.63, 0.76, "#mu_{ref} = "+str(int(cnst.mu_2))+" GeV")
    latexLabel2.DrawLatex(0.63, 0.71, "#mu_{0} = m_{t}")
    if preliminary: latexLabel2.DrawLatex(0.205, 0.92 , "#it{Preliminary}")
    
    c.SaveAs(outdir+'/test_incl.png')
    c.SaveAs(outdir+'/test_incl.pdf')
    c.SaveAs(outdir+'/test_incl.root')
    c.SaveAs(outdir+'/test_incl.C')


    graph.Clear()

    graph.SetPoint(0,cnst.mu_1,ratio_12/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_1, 'nominal')/mass_2))
    graph.SetPointError(0,0,0,err_12_down/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_1, 'nominal')/mass_2),err_12_up/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_1, 'nominal')/mass_2))

    graph.SetPoint(1,cnst.mu_2,1)
    graph.SetPointError(1,0,0,0,0)

    graph.SetPoint(2,cnst.mu_3,ratio_32/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_3, 'nominal')/mass_2))
    graph.SetPointError(2,0,0,err_32_down/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_3, 'nominal')/mass_2),err_32_up/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_3, 'nominal')/mass_2))

    graph.SetPoint(3,cnst.mu_4,ratio_42/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_4, 'nominal')/mass_2))
    graph.SetPointError(3,0,0,err_42_down/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_4, 'nominal')/mass_2),err_42_up/(mtmu2mtmu(mass_2, cnst.mu_2, cnst.mu_4, 'nominal')/mass_2))

    line = TLine(370,1,1030,1)
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
    latexLabel2.DrawLatex(0.59, 0.8, "ABMP16_5_nlo PDF set")

    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/mt_mt.png')
    c.SaveAs(outdir+'/mt_mt.pdf')
    c.SaveAs(outdir+'/mt_mt.root')

    
    gr_add.Clear()

    gr_add.SetPoint(0,cnst.mtmt,cnst.mtmt/mtmt_2)
    gr_add.SetPointError(0,0,cnst.mtmt_err/mtmt_2)

    l = makeAdditionalTheoryPrediction(cnst.mtmt, cnst.mtmt_err, cnst.mtmt_err, mtmu, True)
    scales = l[3]

    gr_band.Clear()
    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],(cnst.mtmt+cnst.mtmt_err)/mtmt_2)
        gr_band.SetPoint(len(scales)+i,scales[len(scales)-i-1],(cnst.mtmt-cnst.mtmt_err)/mtmt_2)

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(m_{t}) / m_{t}^{#mu_{2}}(m_{t})')
    gr_band.SetTitle('m_{t}(m_{t}) measured as a function of #mu')

    gr_band.SetFillStyle(3002)
    gr_band.SetFillColor(rt.kRed+1)
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
    gr_band.SetMinimum(0.88)
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

# "getTotalMassError" calculates the total error in the extracted masses

################################       

def getTotalMassError(mtmu_1, err_1, mtmu_2, err_2, mtmu_3, err_3, mtmu_4, err_4):

    pdf_error_1_up = pdf_error_1_down = 0
    extr_error_1_up = extr_error_1_down = 0
    scale_error_1_up = scale_error_1_down = mtmu_1
    pdf_error_2_up = pdf_error_2_down = 0
    extr_error_2_up = extr_error_2_down = 0
    scale_error_2_up = scale_error_2_down = mtmu_2
    pdf_error_3_up = pdf_error_3_down = 0
    extr_error_3_up = extr_error_3_down = 0
    scale_error_3_up = scale_error_3_down = mtmu_3
    pdf_error_4_up = pdf_error_4_down = 0
    extr_error_4_up = extr_error_4_down = 0
    scale_error_4_up = scale_error_4_down = mtmu_4

    for pdf in range(1,30):
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', pdf , 0 , 0 )

        if mass_and_err_1[0] > mtmu_1 : pdf_error_1_up += (mass_and_err_1[0]-mtmu_1)**2
        else: pdf_error_1_down += (mass_and_err_1[0]-mtmu_1)**2
        if mass_and_err_2[0] > mtmu_2 : pdf_error_2_up += (mass_and_err_2[0]-mtmu_2)**2
        else: pdf_error_2_down += (mass_and_err_2[0]-mtmu_2)**2
        if mass_and_err_3[0] > mtmu_3 : pdf_error_3_up += (mass_and_err_3[0]-mtmu_3)**2
        else: pdf_error_3_down += (mass_and_err_3[0]-mtmu_3)**2
        if mass_and_err_4[0] > mtmu_4 : pdf_error_4_up += (mass_and_err_4[0]-mtmu_4)**2
        else: pdf_error_4_down += (mass_and_err_4[0]-mtmu_4)**2
        
    pdf_error_1_up = pdf_error_1_up**.5
    pdf_error_1_down = pdf_error_1_down**.5
    pdf_error_2_up = pdf_error_2_up**.5
    pdf_error_2_down = pdf_error_2_down**.5
    pdf_error_3_up = pdf_error_3_up**.5
    pdf_error_3_down = pdf_error_3_down**.5
    pdf_error_4_up = pdf_error_4_up**.5
    pdf_error_4_down = pdf_error_4_down**.5

    for extr in range(-len(var.extr_1_up),len(var.extr_1_up)+1) :
        if extr == 0 : continue        
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , extr , 0 )

        if mass_and_err_1[0] > mtmu_1 : extr_error_1_up += (mass_and_err_1[0]-mtmu_1)**2
        else: extr_error_1_down += (mass_and_err_1[0]-mtmu_1)**2
        if mass_and_err_2[0] > mtmu_2 : extr_error_2_up += (mass_and_err_2[0]-mtmu_2)**2
        else: extr_error_2_down += (mass_and_err_2[0]-mtmu_2)**2
        if mass_and_err_3[0] > mtmu_3 : extr_error_3_up += (mass_and_err_3[0]-mtmu_3)**2
        else: extr_error_3_down += (mass_and_err_3[0]-mtmu_3)**2
        if mass_and_err_4[0] > mtmu_4 : extr_error_4_up += (mass_and_err_4[0]-mtmu_4)**2
        else: extr_error_4_down += (mass_and_err_4[0]-mtmu_4)**2

    extr_error_1_up = extr_error_1_up**.5
    extr_error_1_down = extr_error_1_down**.5
    extr_error_2_up = extr_error_2_up**.5
    extr_error_2_down = extr_error_2_down**.5
    extr_error_3_up = extr_error_3_up**.5
    extr_error_3_down = extr_error_3_down**.5
    extr_error_4_up = extr_error_4_up**.5
    extr_error_4_down = extr_error_4_down**.5

    for renscale in 'nominal', 'up','down':
        for facscale in 'nominal', 'up','down':

            if renscale == 'up' and facscale == 'down': continue
            if facscale == 'up' and renscale == 'down': continue
            if facscale == 'nominal' and renscale == 'nominal': continue

            # if facscale_only:
            #     if renscale != 'nominal': continue # only factorization scale
            
            mass_and_err_1 = getMassAndError(1, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_2 = getMassAndError(2, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_3 = getMassAndError(3, renscale, facscale, 0 , 0 , 0 )
            mass_and_err_4 = getMassAndError(4, renscale, facscale, 0 , 0 , 0 )

            # print
            # print 'scales', renscale, facscale, mass_and_err_1[0] - mtmu_1,  mass_and_err_2[0] - mtmu_2, mass_and_err_3[0] - mtmu_3, mass_and_err_4[0] - mtmu_4
            
            if mass_and_err_1[0] > scale_error_1_up : scale_error_1_up = mass_and_err_1[0]
            elif mass_and_err_1[0] < scale_error_1_down : scale_error_1_down = mass_and_err_1[0]
            if mass_and_err_2[0] > scale_error_2_up : scale_error_2_up = mass_and_err_2[0]
            elif mass_and_err_2[0] < scale_error_2_down : scale_error_2_down = mass_and_err_2[0]
            if mass_and_err_3[0] > scale_error_3_up : scale_error_3_up = mass_and_err_3[0]
            elif mass_and_err_3[0] < scale_error_3_down : scale_error_3_down = mass_and_err_3[0]
            if mass_and_err_4[0] > scale_error_4_up : scale_error_4_up = mass_and_err_4[0]
            elif mass_and_err_4[0] < scale_error_4_down : scale_error_4_down = mass_and_err_4[0]
                
    # automatically zero if no scale variations
    scale_error_1_up = scale_error_1_up - mtmu_1
    scale_error_1_down = mtmu_1 - scale_error_1_down
    scale_error_2_up = scale_error_2_up - mtmu_2
    scale_error_2_down = mtmu_2 - scale_error_2_down
    scale_error_3_up = scale_error_3_up - mtmu_3
    scale_error_3_down = mtmu_3 - scale_error_3_down
    scale_error_4_up = scale_error_4_up - mtmu_4
    scale_error_4_down = mtmu_4 - scale_error_4_down
    
    fit_error_1_up = (pdf_error_1_up**2 + extr_error_1_up**2 + err_1**2)**.5
    fit_error_1_down = (pdf_error_1_down**2 + extr_error_1_down**2 + err_1**2)**.5
    fit_error_2_up = (pdf_error_2_up**2 + extr_error_2_up**2 + err_2**2)**.5
    fit_error_2_down = (pdf_error_2_down**2 + extr_error_2_down**2 + err_2**2)**.5
    fit_error_3_up = (pdf_error_3_up**2 + extr_error_3_up**2 + err_3**2)**.5
    fit_error_3_down = (pdf_error_3_down**2 + extr_error_3_down**2 + err_3**2)**.5
    fit_error_4_up = (pdf_error_4_up**2 + extr_error_4_up**2 + err_4**2)**.5
    fit_error_4_down = (pdf_error_4_down**2 + extr_error_4_down**2 + err_4**2)**.5

    tot_error_1_up = (pdf_error_1_up**2 + extr_error_1_up**2 + scale_error_1_up**2 + err_1**2)**.5
    tot_error_1_down = (pdf_error_1_down**2 + extr_error_1_down**2 + scale_error_1_down**2 + err_1**2)**.5
    tot_error_2_up = (pdf_error_2_up**2 + extr_error_2_up**2 + scale_error_2_up**2 + err_2**2)**.5
    tot_error_2_down = (pdf_error_2_down**2 + extr_error_2_down**2 + scale_error_2_down**2 + err_2**2)**.5
    tot_error_3_up = (pdf_error_3_up**2 + extr_error_3_up**2 + scale_error_3_up**2 + err_3**2)**.5
    tot_error_3_down = (pdf_error_3_down**2 + extr_error_3_down**2 + scale_error_3_down**2 + err_3**2)**.5
    tot_error_4_up = (pdf_error_4_up**2 + extr_error_4_up**2 + scale_error_4_up**2 + err_4**2)**.5
    tot_error_4_down = (pdf_error_4_down**2 + extr_error_4_down**2 + scale_error_4_down**2 + err_4**2)**.5


    print
    print 'mt_mu1 =', round(mtmu_1,1), '+/-', round(err_1,1), '(fit) +', round(pdf_error_1_up,1), '-', round(pdf_error_1_down,1), '(pdf) +', round(extr_error_1_up,1), '-', round(extr_error_1_down,1), '(extr) +', round(scale_error_1_up,1), '-', round(scale_error_1_down,1), '(scale) =', round(mtmu_1,1), '+', round(tot_error_1_up,1), '-', round(tot_error_1_down,1), '(tot)' 
    print 'mt_mu1 =', round(mtmu_2,1), '+/-', round(err_2,1), '(fit) +', round(pdf_error_2_up,1), '-', round(pdf_error_2_down,1), '(pdf) +', round(extr_error_2_up,1), '-', round(extr_error_2_down,1), '(extr) +', round(scale_error_2_up,1), '-', round(scale_error_2_down,1), '(scale) =', round(mtmu_2,1), '+', round(tot_error_2_up,1), '-', round(tot_error_2_down,1), '(tot)' 
    print 'mt_mu3 =', round(mtmu_3,1), '+/-', round(err_3,1), '(fit) +', round(pdf_error_3_up,1), '-', round(pdf_error_3_down,1), '(pdf) +', round(extr_error_3_up,1), '-', round(extr_error_3_down,1), '(extr) +', round(scale_error_3_up,1), '-', round(scale_error_3_down,1), '(scale) =', round(mtmu_3,1), '+', round(tot_error_3_up,1), '-', round(tot_error_3_down,1), '(tot)' 
    print 'mt_mu4 =', round(mtmu_4,1), '+/-', round(err_4,1), '(fit) +', round(pdf_error_4_up,1), '-', round(pdf_error_4_down,1), '(pdf) +', round(extr_error_4_up,1), '-', round(extr_error_4_down,1), '(extr) +', round(scale_error_4_up,1), '-', round(scale_error_4_down,1), '(scale) =', round(mtmu_4,1), '+', round(tot_error_4_up,1), '-', round(tot_error_4_down,1), '(tot)' 
    print
    
    return [fit_error_1_up, fit_error_1_down, fit_error_2_up, fit_error_2_down, fit_error_3_up, fit_error_3_down, fit_error_4_up, fit_error_4_down,
            tot_error_1_up, tot_error_1_down, tot_error_2_up, tot_error_2_down, tot_error_3_up, tot_error_3_down, tot_error_4_up, tot_error_4_down]

################################

# "makeMassPlots" produces the absolute mass plots for the running of mt

################################

def makeMassPlots(mtmu_1, err_1, mtmt_1, mtmu_2, err_2, mtmt_2, mtmu_3, err_3, mtmt_3, mtmu_4, err_4, mtmt_4) :

    errors = getTotalMassError(mtmu_1, err_1, mtmu_2, err_2, mtmu_3, err_3, mtmu_4, err_4)

    graph=TGraphAsymmErrors(3)
    graph.SetPoint(0,cnst.mu_1,mtmu_1)
    graph.SetPointError(0,0,0,err_1,err_1)
    graph.SetPoint(1,cnst.mu_2,mtmu_2)
    graph.SetPointError(1,0,0,err_2,err_2)
    graph.SetPoint(2,cnst.mu_3,mtmu_3)
    graph.SetPointError(2,0,0,err_3,err_3)
    graph.SetPoint(3,cnst.mu_4,mtmu_4)
    graph.SetPointError(3,0,0,err_4,err_4)

    gr_add=TGraphErrors()
    gr_add.SetPoint(0,cnst.mtmt,cnst.mtmt)
    gr_add.SetPointError(0,0,cnst.mtmt_err)

    
    l = makeAdditionalTheoryPrediction(cnst.mtmt, cnst.mtmt_err, cnst.mtmt_err, mtmu_2, False)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]
    
    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1])

    gr_band.SetFillStyle(3002)
    gr_band.SetFillColor(rt.kRed+1)

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(#mu) [GeV]')
    gr_band.SetTitle('')

    gr_add.SetMarkerStyle(21)
    gr_add.SetMarkerColor(rt.kBlue)
    gr_add.SetLineColor(rt.kBlue)

    gr_band.GetYaxis().SetRangeUser(140,170)

    graph.SetMarkerStyle(8)
    graph.SetTitle('')
    
    # leg = TLegend(.15,.23,.7,.35)
    leg = TLegend(.15,.15,.77,.3)
    leg.SetBorderSize(0)
    leg.AddEntry(graph,'NLO extraction from differential #sigma_{t#bar{t}}','pe')
    leg.AddEntry(gr_add,'NLO extraction from inclusive #sigma_{t#bar{t}} (same data)','pe')
    leg.AddEntry(gr_band,'Evolved uncertainty, one loop RGE (5 flavours)','f')
    
    c = TCanvas()
    gr_band.SetMinimum(125)
    gr_band.Draw('af')
    gr_add.Draw('p same')
    graph.Draw('p same')
    leg.Draw('same')

    #fromhere
    latexLabel1 = TLatex()
    latexLabel1.SetTextSize(0.06)
    latexLabel1.SetNDC()
    
    latexLabel2 = TLatex()
    latexLabel2.SetTextSize(0.04)
    latexLabel2.SetTextFont(42)
    latexLabel2.SetNDC()

    latexLabel1.DrawLatex(0.11, 0.92, "CMS")
    latexLabel2.DrawLatex(0.70, 0.92, "35.9 fb^{-1} (13 TeV)")
    latexLabel2.DrawLatex(0.59, 0.78, "ABMP16_5_nlo PDF set")
    # latexLabel2.DrawLatex(0.59, 0.73, "#mu_{ref} = "+str(round(cnst.mu_2,1))+" GeV")

    
    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/test_mtmu_abs.png')
    c.SaveAs(outdir+'/test_mtmu_abs.pdf')
    c.SaveAs(outdir+'/test_mtmu_abs.root')

    graph.SetPointError(0,0,0,errors[1],errors[0])
    graph.SetPointError(1,0,0,errors[3],errors[2])
    graph.SetPointError(2,0,0,errors[5],errors[4])
    graph.SetPointError(3,0,0,errors[7],errors[6])

    c.Update()
    c.SaveAs(outdir+'/test_mtmu_abs_uncert.png')
    c.SaveAs(outdir+'/test_mtmu_abs_uncert.pdf')
    c.SaveAs(outdir+'/test_mtmu_abs_uncert.root')
    
    l = makeAdditionalTheoryPrediction(cnst.mtmt, (cnst.mtmt_err**2 + cnst.mtmt_scale_up**2)**.5, (cnst.mtmt_err**2 + cnst.mtmt_scale_down**2)**.5, mtmu_2, False)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]
    
    gr_band_fit = gr_band.Clone()
    gr_band_fit.SetFillColor(rt.kBlue-1)
    gr_band.Clear()

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1])

    # gr_band.SetFillStyle(3002)
    # gr_band.SetFillColor(rt.kRed+1)

    graph_scale = graph.Clone()
    graph_scale.SetPointError(0,0,0,errors[9],errors[8])
    graph_scale.SetPointError(1,0,0,errors[11],errors[10])
    graph_scale.SetPointError(2,0,0,errors[13],errors[12])
    graph_scale.SetPointError(3,0,0,errors[15],errors[14])
    graph.SetLineWidth(2)
    gr_add.SetLineWidth(2)
    
    gr_add_scale=TGraphAsymmErrors()
    gr_add_scale.SetPoint(0,cnst.mtmt,cnst.mtmt)
    gr_add_scale.SetPointError(0,0,0,(cnst.mtmt_err**2+cnst.mtmt_scale_down**2)**.5,(cnst.mtmt_err**2+cnst.mtmt_scale_up**2)**.5)
    gr_add_scale.SetMarkerStyle(7)
    gr_add_scale.SetMarkerColor(rt.kBlue)
    gr_add_scale.SetLineColor(rt.kBlue)

    
    gr_band.SetMinimum(110)
    gr_band.SetMaximum(175)
    graph_scale.Draw('p same')
    gr_add_scale.Draw('p same')
    gr_band_fit.Draw('f same')
    c.Update()
    c.SaveAs(outdir+'/test_mtmu_abs_scale.png')
    c.SaveAs(outdir+'/test_mtmu_abs_scale.pdf')
    c.SaveAs(outdir+'/test_mtmu_abs_scale.root')
    
    

    graph.Clear()
    graph=TGraphErrors(3)
    graph.SetPoint(0,cnst.mu_1,mtmt_1)
    graph.SetPointError(0,0,err_1)
    graph.SetPoint(1,cnst.mu_2,mtmt_2)
    graph.SetPointError(1,0,err_2)
    graph.SetPoint(2,cnst.mu_3,mtmt_3)
    graph.SetPointError(2,0,err_3)
    graph.SetPoint(3,cnst.mu_4,mtmt_4)
    graph.SetPointError(3,0,err_4)
    
    gr_band.Clear()
    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],cnst.mtmt+cnst.mtmt_err)
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],cnst.mtmt-cnst.mtmt_err)

    gr_band.SetFillStyle(3002)
    gr_band.SetFillColor(rt.kRed+1)

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(m_{t}) [GeV]')
    gr_band.SetTitle('m_{t}(m_{t}) measured as a function of #mu')

    gr_add.SetMarkerStyle(21)
    gr_add.SetMarkerColor(rt.kBlue)
    gr_add.SetLineColor(rt.kBlue)

    gr_band.GetYaxis().SetRangeUser(157,170)

    graph.SetMarkerStyle(8)
    
    leg.Clear()
    leg = TLegend(.15,.23,.7,.35)
    leg.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}} (exp. only)','pe')
    leg.AddEntry(gr_add,'Hathor @NLO from incl. #sigma_{t#bar{t}} (same data)')
    leg.AddEntry(gr_band,'experimental + extrapolation + PDF uncertainties','f')
    
    c.Clear()
    gr_band.SetMinimum(142)
    gr_band.SetMaximum(170)
    gr_band.Draw('af')
    gr_add.Draw('p same')
    graph.Draw('p same')
    leg.Draw('same')

    c.SaveAs(outdir+'/test_mtmt_abs.png')
    c.SaveAs(outdir+'/test_mtmt_abs.pdf')
    c.SaveAs(outdir+'/test_mtmt_abs.root')
  
    return


################################

# "makeChi2Test" calculates the significance of the obseved running
#  using a chi2 test

################################


def makeChi2Test(mass2, ratio12, ratio32, ratio42, err12_up, err12_down, err32_up, err32_down, err42_up, err42_down):

    th_ratio12 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_1, 'nominal')/mass2
    th_ratio32 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_3, 'nominal')/mass2
    th_ratio42 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_4, 'nominal')/mass2

    print 
    print 'theory:', round(th_ratio12,3), round(th_ratio32,3), round(th_ratio42,3)

    if ratio12 > th_ratio12: err12 = err12_down
    else: err12 = err12_up

    if ratio32 > th_ratio32: err32 = err32_down
    else: err32 = err32_up

    if ratio42 > th_ratio42: err42 = err42_down
    else: err42 = err42_up

    chi2_no_run = ((ratio12-1)/err12)**2 + ((ratio32-1)/err32)**2 + ((ratio42-1)/err42)**2
    chi2_run = ((ratio12-th_ratio12)/err12)**2 + ((ratio32-th_ratio32)/err32)**2 + ((ratio42-th_ratio42)/err42)**2

    print 
    print 'chi2: running hp', round(chi2_run,2)
    print 'chi2: no running hp', round(chi2_no_run,2)
    print
    print 'significance =', round((chi2_no_run-chi2_run)**.5,2)
    print
    
    
    return

################################

# "getChi2Uncorrelated" performs a transformation to an orthogonal basis
#  before calculating the chisquare, in order to take correlations into account.

################################


def getChi2Uncorrelated(ratio12, th1, err12, ratio32, th3, err32, ratio42, th4, err42):
    diff_1_2 = (ratio12 - th1)
    diff_3_2 = (ratio32 - th3)
    diff_4_2 = (ratio42 - th4)

    V = [[err12*err12, err12*err32*rcorr_12, err12*err42*rcorr_13],
         [err12*err32*rcorr_12, err32*err32, err32*err42*rcorr_23],
         [err12*err42*rcorr_13, err32*err42*rcorr_23, err42*err42],]

    E = np.linalg.eigh(V)
    values = E[0]
    vectors = E[1]
    
    x = np.array([diff_1_2, diff_3_2, diff_4_2])
    y1 = np.inner(np.array(vectors[0]),x)
    y2 = np.inner(np.array(vectors[1]),x)
    y3 = np.inner(np.array(vectors[2]),x)
    
    chi2 = (y1**2)/values[0] + (y2**2)/values[1] + (y3**2)/values[2]
    return chi2

    
################################

# "makeChi2Significance" calculates the significance of the obseved running
#  by fitting the data to the theory prediction

################################


def makeChi2Significance(mass2, ratio12, ratio32, ratio42, err12, err32, err42):

    th_ratio12 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_1, 'nominal')/mass2
    th_ratio32 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_3, 'nominal')/mass2
    th_ratio42 = mtmu2mtmu(mass2, cnst.mu_2, cnst.mu_4, 'nominal')/mass2

    global ntoys
    if ntoys < 10000 :
        if ntoys < 1000: ntoys = 1000
        print '\nWARNING: running only', ntoys ,'toys for estimate of correlations.'
        print '         It is recommended to run at least 10.000'
    estimateMassCorrelations()
        
    graph = TGraph()
    
    for x in range(0,31):
        x = x/10.
        th1 = x*(th_ratio12-1)+1
        th3 = x*(th_ratio32-1)+1
        th4 = x*(th_ratio42-1)+1
        chi2 = getChi2Uncorrelated(ratio12, th1, err12, ratio32, th3, err32, ratio42, th4, err42)
        
        graph.SetPoint(int(x*10),x,chi2)

    funct = TF1('funct','pol2(0)',0,3)
    graph.Fit(funct,'q','',0,3)

    xmin = funct.GetMinimumX()
    ymin = funct.GetMinimum(0,3)
    err = xmin - funct.GetX(ymin+1,0,xmin)
   
    c = TCanvas()
    graph.Draw('ap')
    graph.SetMarkerStyle(8)
    graph.SetMinimum(0)

    graph.GetXaxis().SetTitle('x')
    graph.GetYaxis().SetTitle('#chi^{2}')
    
    outdir = 'plots_chi2'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.Print(outdir+'/chi2_significance_nominal.png')


    err_pdf_up = 0
    err_pdf_down = 0
    
    for pdf in range(1,30):
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', pdf , 0 , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', pdf , 0 , 0 )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]
        ratio_4_2 = ratios_and_errs[2]
        err_ratio_1_2 = ratios_and_errs[3]
        err_ratio_3_2 = ratios_and_errs[4]
        err_ratio_4_2 = ratios_and_errs[5]

        graph.Clear()
        
        for x in range(0,31):
            x = x/10.
            th1 = x*(th_ratio12-1)+1
            th3 = x*(th_ratio32-1)+1
            th4 = x*(th_ratio42-1)+1
            chi2 = getChi2Uncorrelated(ratio_1_2, th1, err_ratio_1_2, ratio_3_2, th3, err_ratio_3_2, ratio_4_2, th4, err_ratio_4_2)
            graph.SetPoint(int(x*10),x,chi2)

        funct = TF1('funct','pol2(0)',0,3)
        graph.Fit(funct,'q','',0,3)

        c.Clear()
        graph.Draw('ap')
        c.Print(outdir+'/chi2_significance_PDF'+str(pdf)+'.png')
        
        xmin_PDF = funct.GetMinimumX()
        err_pdf = (xmin-xmin_PDF)**2
        
        if xmin_PDF > xmin : err_pdf_up += err_pdf
        else : err_pdf_down += err_pdf

    err_pdf_up = err_pdf_up**.5
    err_pdf_down = err_pdf_up**.5

    err_extr_up = 0
    err_extr_down = 0

    for extr in range(-len(var.extr_1_up),len(var.extr_1_up)+1) :
        if extr == 0 : continue

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , extr , 0 )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , extr , 0 )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]
        ratio_4_2 = ratios_and_errs[2]
        err_ratio_1_2 = ratios_and_errs[3]
        err_ratio_3_2 = ratios_and_errs[4]
        err_ratio_4_2 = ratios_and_errs[5]


        graph.Clear()
        
        for x in range(0,31):
            x = x/10.
            th1 = x*(th_ratio12-1)+1
            th3 = x*(th_ratio32-1)+1
            th4 = x*(th_ratio42-1)+1
            chi2 = getChi2Uncorrelated(ratio_1_2, th1, err_ratio_1_2, ratio_3_2, th3, err_ratio_3_2, ratio_4_2, th4, err_ratio_4_2)
            graph.SetPoint(int(x*10),x,chi2)

        funct = TF1('funct','pol2(0)',0,3)
        graph.Fit(funct,'q','',0,3)

        c.Clear()
        graph.Draw('ap')
        c.Print(outdir+'/chi2_significance_EXTR'+str(extr)+'.png')
        
        xmin_EXTR = funct.GetMinimumX()
        err_extr = (xmin-xmin_EXTR)**2
        
        if xmin_EXTR > xmin : err_extr_up += err_extr
        else : err_extr_down += err_extr

    err_extr_up = err_extr_up**.5
    err_extr_down = err_extr_up**.5

    err_scale_up = 0
    err_scale_down = 0

    if do_scale_variations:
        for facscale in ['up','down']:
            mass_and_err_1 = getMassAndError(1, 'nominal', facscale, 0 , 0 , 0 )
            mass_and_err_2 = getMassAndError(2, 'nominal', facscale, 0 , 0 , 0 )
            mass_and_err_3 = getMassAndError(3, 'nominal', facscale, 0 , 0 , 0 )
            mass_and_err_4 = getMassAndError(4, 'nominal', facscale, 0 , 0 , 0 )

            ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                        mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])
            
            ratio_1_2 = ratios_and_errs[0]
            ratio_3_2 = ratios_and_errs[1]
            ratio_4_2 = ratios_and_errs[2]
            err_ratio_1_2 = ratios_and_errs[3]
            err_ratio_3_2 = ratios_and_errs[4]
            err_ratio_4_2 = ratios_and_errs[5]

            graph.Clear()
        
            for x in range(0,31):
                x = x/10.
                th1 = x*(th_ratio12-1)+1
                th3 = x*(th_ratio32-1)+1
                th4 = x*(th_ratio42-1)+1
                chi2 = getChi2Uncorrelated(ratio_1_2, th1, err_ratio_1_2, ratio_3_2, th3, err_ratio_3_2, ratio_4_2, th4, err_ratio_4_2)
                graph.SetPoint(int(x*10),x,chi2)

            funct = TF1('funct','pol2(0)',0,3)
            graph.Fit(funct,'q','',0,3)

            c.Clear()
            graph.Draw('ap')
            c.Print(outdir+'/chi2_significance_facscale_'+facscale+'.png')
            
            xmin_SCALE = funct.GetMinimumX()
            err_scale = (xmin-xmin_SCALE)**2
        
            if xmin_SCALE > xmin : err_scale_up += err_scale
            else : err_scale_down += err_scale

        err_scale_up = err_scale_up**.5
        err_scale_down = err_scale_up**.5
        
        err_up = (err**2 + err_pdf_up**2 + err_extr_up**2 + err_scale_up**2)**.5
        err_down = (err**2 + err_pdf_down**2 + err_extr_down**2 + err_scale_down**2)**.5

    else:
        err_up = (err**2 + err_pdf_up**2 + err_extr_up**2)**.5
        err_down = (err**2 + err_pdf_down**2 + err_extr_down**2)**.5



        
    print
    if do_scale_variations:
        print 'xmin =', round(xmin,2), '+/-' ,round(err,2), '(exp)', '+' ,round(err_pdf_up,2), '-', round(err_pdf_down,2), '(PDF)', '+' ,round(err_extr_up,2), '-', round(err_extr_down,2), '(extr)' , '+' ,round(err_scale_up,2), '-', round(err_scale_down,2), '(scale)'
    else : print 'xmin =', round(xmin,2), '+/-' ,round(err,2), '(exp)', '+' ,round(err_pdf_up,2), '-', round(err_pdf_down,2), '(PDF)', '+' ,round(err_extr_up,2), '-', round(err_extr_down,2), '(extr)'
    print
    print 'xmin =', round(xmin,2), '+' ,round(err_up,2), '-', round(err_down,2), '(tot)'
    print
    print 'significance wrt no running =', round(xmin/err_down,2)
    print 'significance wrt RGE =', round((xmin-1)/err_down,2)
    print




    
    return




################################

# "estimateSystContributions" estimates the contribution to the total uncertainty
#  of each individual systematics source

################################


def estimateSystContributions(central_ratio_1_2, central_ratio_3_2, central_ratio_4_2):
    
    print
    print 'contribs\n'

    for contrib in range(1,len(var.contribs_1)+1) :

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , contrib )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , contrib )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , contrib )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , 0 , contrib )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2_up = ratios_and_errs[0]
        ratio_3_2_up = ratios_and_errs[1]
        ratio_4_2_up = ratios_and_errs[2]

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , 0 , int(-1*contrib) )

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        ratio_1_2_down = ratios_and_errs[0]
        ratio_3_2_down = ratios_and_errs[1]
        ratio_4_2_down = ratios_and_errs[2]

        contr_1_2_up = (ratio_1_2_up/central_ratio_1_2-1)*100
        contr_1_2_down = (ratio_1_2_down/central_ratio_1_2-1)*100
        contr_3_2_up = (ratio_3_2_up/central_ratio_3_2-1)*100
        contr_3_2_down = (ratio_3_2_down/central_ratio_3_2-1)*100
        contr_4_2_up = (ratio_4_2_up/central_ratio_4_2-1)*100
        contr_4_2_down = (ratio_4_2_down/central_ratio_4_2-1)*100

        
        if (contr_1_2_up*contr_1_2_down) < 0:
            contr_1_2 = abs(contr_1_2_up-contr_1_2_down)/2.
        else : contr_1_2 = max(abs(contr_1_2_up),abs(contr_1_2_down))

        if (contr_3_2_up*contr_3_2_down) < 0:
            contr_3_2 = abs(contr_3_2_up-contr_3_2_down)/2.
        else : contr_3_2 = max(abs(contr_3_2_up),abs(contr_3_2_down))

        if (contr_4_2_up*contr_4_2_down) < 0:
            contr_4_2 = abs(contr_4_2_up-contr_4_2_down)/2.
        else : contr_4_2 = max(abs(contr_4_2_up),abs(contr_4_2_down))

        
        if contr_1_2 > 0.2 or contr_3_2 > 0.3 or contr_4_2 > 0.5:
            print var.syst_names[contrib-1], round(contr_1_2,2), round(contr_3_2,2), round(contr_4_2,2)

    print
    
    return



################################

# "throwToyCrossSections" throws toy cross sections taking the correlations into account
#  the results are used to estimate the correlation between the fitted masses

################################

def throwToyCrossSections(r):

    global xsec_1, xsec_2, xsec_3, xsec_4

    err_1 = xsec_1*(var.err_xsec_1_up+var.err_xsec_1_down)/2./100.
    err_2 = xsec_2*(var.err_xsec_2_up+var.err_xsec_2_down)/2./100.
    err_3 = xsec_3*(var.err_xsec_3_up+var.err_xsec_3_down)/2./100.
    err_4 = xsec_4*(var.err_xsec_4_up+var.err_xsec_4_down)/2./100.

    V = [[err_1*err_1, var.corr_1_2*err_1*err_2, var.corr_1_3*err_1*err_3, var.corr_4_1*err_1*err_4],
         [var.corr_1_2*err_1*err_2, err_2*err_2, var.corr_3_2*err_2*err_3, var.corr_4_2*err_2*err_4],
         [var.corr_1_3*err_1*err_3, var.corr_3_2*err_2*err_3, err_3*err_3, var.corr_4_3*err_3*err_4],
         [var.corr_4_1*err_1*err_4, var.corr_4_2*err_2*err_4, var.corr_4_3*err_3*err_4, err_4*err_4]]



    L = np.linalg.cholesky(V) # Cholesky decomposition 

    mu = [xsec_1, xsec_2, xsec_3, xsec_4]
    z = [r.Gaus(0,1),r.Gaus(0,1),r.Gaus(0,1),r.Gaus(0,1)]

    y = L.dot(z)
    y += mu

    xsec_1 = y[0]
    xsec_2 = y[1]
    xsec_3 = y[2]
    xsec_4 = y[3]
        
    return



################################

# "estimateMassCorrelations" estimates the correlation between the fitted masses
#  to run, set "ntoys" to some positive numbers (at least 10k)
#  N.B. it is absolutely safe to use the correlaton between the cross section (i.e. ntoys=0)

################################


def estimateMassCorrelations():

    global doingToys
    global xsec_1, xsec_2, xsec_3, xsec_4
    
    orig_xsec_1 = xsec_1
    orig_xsec_2 = xsec_2
    orig_xsec_3 = xsec_3
    orig_xsec_4 = xsec_4
    
    doingToys = True
    
    r=TRandom3()

    m12 = TH2D('m12','m12',100,150,160,100,130,170)
    m32 = TH2D('m32','m32',100,120,180,100,130,170)
    m42 = TH2D('m42','m42',100,110,170,100,130,170)

    xs12 = TH2D('xs12','xs12',100,xsec_1*(1-5*var.err_xsec_1_down/100.),xsec_1*(1+5*var.err_xsec_1_up/100.),
                100,xsec_2*(1-5*var.err_xsec_2_down/100.),xsec_2*(1-5*var.err_xsec_2_down/100.))
    xs13 = TH2D('xs13','xs13',100,xsec_1*(1-5*var.err_xsec_1_down/100.),xsec_1*(1+5*var.err_xsec_1_up/100.),
                100,xsec_3*(1-5*var.err_xsec_3_down/100.),xsec_3*(1-5*var.err_xsec_3_down/100.))
    xs23 = TH2D('xs23','xs23',100,xsec_2*(1-5*var.err_xsec_2_down/100.),xsec_2*(1+5*var.err_xsec_2_up/100.),
                100,xsec_3*(1-5*var.err_xsec_3_down/100.),xsec_3*(1-5*var.err_xsec_3_down/100.))
    xs24 = TH2D('xs24','xs24',100,xsec_2*(1-5*var.err_xsec_2_down/100.),xsec_2*(1+5*var.err_xsec_2_up/100.),
                100,xsec_4*(1-5*var.err_xsec_4_down/100.),xsec_4*(1-5*var.err_xsec_4_down/100.))
    xs14 = TH2D('xs14','xs14',100,xsec_1*(1-5*var.err_xsec_1_down/100.),xsec_1*(1+5*var.err_xsec_1_up/100.),
                100,xsec_4*(1-5*var.err_xsec_4_down/100.),xsec_4*(1-5*var.err_xsec_4_down/100.))
    xs34 = TH2D('xs34','xs34',100,xsec_3*(1-5*var.err_xsec_3_down/100.),xsec_3*(1+5*var.err_xsec_3_up/100.),
                100,xsec_4*(1-5*var.err_xsec_4_down/100.),xsec_4*(1-5*var.err_xsec_4_down/100.))


    
    r_corr_12 = TH2D('r_corr_12','r_corr_12',100,0.9,1.15,100,0.75,1.15)
    r_corr_13 = TH2D('r_corr_13','r_corr_13',100,0.9,1.15,100,0.75,1.15)
    r_corr_23 = TH2D('r_corr_23','r_corr_23',100,0.9,1.15,100,0.75,1.15)
    
    print '\nexecuting', ntoys, 'toys\n'

    for i in range(1,ntoys+1):
        if i%1000==0 or i==1: print 'toy n.', i
        throwToyCrossSections(r)

        xs12.Fill(xsec_1,xsec_2)
        xs13.Fill(xsec_1,xsec_3)
        xs23.Fill(xsec_2,xsec_3)
        xs14.Fill(xsec_1,xsec_4)
        xs24.Fill(xsec_2,xsec_4)
        xs34.Fill(xsec_3,xsec_4)
        
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , 0)
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , 0)
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , 0)
        mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , 0 , 0)
        
        m12.Fill(mass_and_err_1[0],mass_and_err_2[0])
        m32.Fill(mass_and_err_3[0],mass_and_err_2[0])
        m42.Fill(mass_and_err_4[0],mass_and_err_2[0])

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

        r_corr_12.Fill(ratios_and_errs[0],ratios_and_errs[1])
        r_corr_13.Fill(ratios_and_errs[0],ratios_and_errs[2])
        r_corr_23.Fill(ratios_and_errs[1],ratios_and_errs[2])
        
        #extremely important
        xsec_1 = orig_xsec_1
        xsec_2 = orig_xsec_2
        xsec_3 = orig_xsec_3
        xsec_4 = orig_xsec_4
    
    doingToys = False

    if replace_corr:
        var.corr_1_2 = m12.GetCorrelationFactor()
        var.corr_3_2 = m32.GetCorrelationFactor()
        var.corr_4_2 = m42.GetCorrelationFactor()
    
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
    m42.SetTitle('correlation = '+str(round(m42.GetCorrelationFactor(),2)))
    m42.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_4))+' GeV) [GeV]')
    m42.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_2))+' GeV) [GeV]')
    m42.DrawNormalized('colz')
    c.SaveAs(outdir+'/mass_corr_4_2.png')
    c.SaveAs(outdir+'/mass_corr_4_2.pdf')
    c.SaveAs(outdir+'/mass_corr_4_2.root')
    c.Clear()
    r_corr_12.SetTitle('correlation = '+str(round(r_corr_12.GetCorrelationFactor(),2)))
    r_corr_12.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_1))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_12.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_3))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_12.DrawNormalized('colz')
    c.SaveAs(outdir+'/ratio_corr_12.png')
    c.SaveAs(outdir+'/ratio_corr_12.pdf')
    c.SaveAs(outdir+'/ratio_corr_12.root')
    r_corr_13.SetTitle('correlation = '+str(round(r_corr_13.GetCorrelationFactor(),2)))
    r_corr_13.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_1))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_13.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_4))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_13.DrawNormalized('colz')
    c.SaveAs(outdir+'/ratio_corr_13.png')
    c.SaveAs(outdir+'/ratio_corr_13.pdf')
    c.SaveAs(outdir+'/ratio_corr_13.root')
    r_corr_23.SetTitle('correlation = '+str(round(r_corr_23.GetCorrelationFactor(),2)))
    r_corr_23.GetXaxis().SetTitle('mt ('+str(int(cnst.mu_3))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_23.GetYaxis().SetTitle('mt ('+str(int(cnst.mu_4))+' GeV) / mt ('+str(int(cnst.mu_2))+' GeV)')
    r_corr_23.DrawNormalized('colz')
    c.SaveAs(outdir+'/ratio_corr_23.png')
    c.SaveAs(outdir+'/ratio_corr_23.pdf')
    c.SaveAs(outdir+'/ratio_corr_23.root')


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
    xs14.SetTitle('correlation = '+str(round(xs14.GetCorrelationFactor(),2)))
    xs14.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{1}) [pb]')
    xs14.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{4}) [pb]')
    xs14.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_14.png')
    c.SaveAs(outdir+'/xsec_corr_14.pdf')
    c.SaveAs(outdir+'/xsec_corr_14.root')
    xs23.SetTitle('correlation = '+str(round(xs23.GetCorrelationFactor(),2)))
    xs23.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{2}) [pb]')
    xs23.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{3}) [pb]')
    xs23.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_23.png')
    c.SaveAs(outdir+'/xsec_corr_23.pdf')
    c.SaveAs(outdir+'/xsec_corr_23.root')
    xs24.SetTitle('correlation = '+str(round(xs24.GetCorrelationFactor(),2)))
    xs24.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{2}) [pb]')
    xs24.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{4}) [pb]')
    xs24.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_24.png')
    c.SaveAs(outdir+'/xsec_corr_24.pdf')
    c.SaveAs(outdir+'/xsec_corr_24.root')
    xs34.SetTitle('correlation = '+str(round(xs34.GetCorrelationFactor(),2)))
    xs34.GetXaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{3}) [pb]')
    xs34.GetYaxis().SetTitle('#sigma_{r#bar{t}} (#mu_{4}) [pb]')
    xs34.DrawNormalized('colz')
    c.SaveAs(outdir+'/xsec_corr_34.png')
    c.SaveAs(outdir+'/xsec_corr_34.pdf')
    c.SaveAs(outdir+'/xsec_corr_34.root')

    global rcorr_12, rcorr_13, rcorr_23
    rcorr_12 = r_corr_12.GetCorrelationFactor()
    rcorr_13 = r_corr_13.GetCorrelationFactor()
    rcorr_23 = r_corr_23.GetCorrelationFactor()
    
    print
    if replace_corr:
        print 'new corr_1_2 =', round(var.corr_1_2,3)
        print 'new corr_3_2 =', round(var.corr_3_2,3)
        print 'new corr_4_2 =', round(var.corr_4_2,3)
    else:    
        print 'estimated corr_1_2 =', round(m12.GetCorrelationFactor(),3)
        print 'estimated corr_3_2 =', round(m32.GetCorrelationFactor(),3)
        print 'estimated corr_4_2 =', round(m42.GetCorrelationFactor(),3)
        print
        
    print 'estimated corr ratios 1_2  =', round(r_corr_12.GetCorrelationFactor(),3)
    print 'estimated corr ratios 1_3  =', round(r_corr_13.GetCorrelationFactor(),3)
    print 'estimated corr ratios 2_3  =', round(r_corr_23.GetCorrelationFactor(),3)
    print
    
    return


def getTotalError (ratios_and_errs, pdf_errors, extr_errors, scale_errors):

    ratio_1_2 = ratios_and_errs[0]
    ratio_3_2 = ratios_and_errs[1]
    ratio_4_2 = ratios_and_errs[2]
    err_ratio_1_2 = ratios_and_errs[3]
    err_ratio_3_2 = ratios_and_errs[4]
    err_ratio_4_2 = ratios_and_errs[5]

    if not newscales: # FIXME #
        err_pdf_1_2_up = pdf_errors[0]
        err_pdf_1_2_down = pdf_errors[1]
        err_pdf_3_2_up = pdf_errors[2]
        err_pdf_3_2_down = pdf_errors[3]
        err_pdf_4_2_up = pdf_errors[4]
        err_pdf_4_2_down = pdf_errors[5]

    err_extr_1_2_up = extr_errors[0]
    err_extr_1_2_down = extr_errors[1]
    err_extr_3_2_up = extr_errors[2]
    err_extr_3_2_down = extr_errors[3]
    err_extr_4_2_up = extr_errors[4]
    err_extr_4_2_down = extr_errors[5]

    if not newscales: # FIXME #
        err_1_2_up = (err_ratio_1_2**2 + err_pdf_1_2_up**2 + err_extr_1_2_up **2)**.5
        err_1_2_down = (err_ratio_1_2**2 + err_pdf_1_2_down**2 + err_extr_1_2_down **2)**.5
        err_3_2_up = (err_ratio_3_2**2 + err_pdf_3_2_up**2 + err_extr_3_2_up **2)**.5
        err_3_2_down = (err_ratio_3_2**2 + err_pdf_3_2_down**2 + err_extr_3_2_down **2)**.5
        err_4_2_up = (err_ratio_4_2**2 + err_pdf_4_2_up**2 + err_extr_4_2_up **2)**.5
        err_4_2_down = (err_ratio_4_2**2 + err_pdf_4_2_down**2 + err_extr_4_2_down **2)**.5
    else:
        err_1_2_up = (err_ratio_1_2**2 + err_extr_1_2_up **2)**.5
        err_1_2_down = (err_ratio_1_2**2 + err_extr_1_2_down **2)**.5
        err_3_2_up = (err_ratio_3_2**2 + err_extr_3_2_up **2)**.5
        err_3_2_down = (err_ratio_3_2**2 + err_extr_3_2_down **2)**.5
        err_4_2_up = (err_ratio_4_2**2 + err_extr_4_2_up **2)**.5
        err_4_2_down = (err_ratio_4_2**2 + err_extr_4_2_down **2)**.5

        
    if do_scale_variations:
        err_scale_1_2_up = scale_errors[0]
        err_scale_1_2_down = scale_errors[1]
        err_scale_3_2_up = scale_errors[2]
        err_scale_3_2_down = scale_errors[3]
        err_scale_4_2_up = scale_errors[4]
        err_scale_4_2_down = scale_errors[5]

        err_1_2_up = (err_1_2_up**2 + err_scale_1_2_up**2)**.5
        err_1_2_down = (err_1_2_down**2 + err_scale_1_2_down**2)**.5
        err_3_2_up = (err_3_2_up**2 + err_scale_3_2_up**2)**.5
        err_3_2_down = (err_3_2_down**2 + err_scale_3_2_down**2)**.5
        err_4_2_up = (err_4_2_up**2 + err_scale_4_2_up**2)**.5
        err_4_2_down = (err_4_2_down**2 + err_scale_4_2_down**2)**.5


    print '\n'
    print 'uncertainties ratio_1_2:\n'
    print 'experimental =', round(err_ratio_1_2,3), round(err_ratio_1_2/ratio_1_2*100.,2), '%'
    if not newscales: # FIXME #
        print 'PDFs up =', round(err_pdf_1_2_up,3), round(err_pdf_1_2_up/ratio_1_2*100.,2), '%'
        print 'PDFs down =', round(err_pdf_1_2_down,3), round(err_pdf_1_2_down/ratio_1_2*100.,2), '%'
    print 'extr up =', round(err_extr_1_2_up,3), round(err_extr_1_2_up/ratio_1_2*100.,2), '%'
    print 'extr down =', round(err_extr_1_2_down,3), round(err_extr_1_2_down/ratio_1_2*100.,2), '%'
    if do_scale_variations:
        print 'scale up =', round(err_scale_1_2_up,3), round(err_scale_1_2_up/ratio_1_2*100.,2), '%'
        print 'scale down =', round(err_scale_1_2_down,3), round(err_scale_1_2_down/ratio_1_2*100.,2), '%'
    print 'total =', round(.5*(err_1_2_up+err_1_2_down)/ratio_1_2*100.,2), '%'

    print '\n'
    print 'uncertainties ratio_3_2:\n'
    print 'experimental =', round(err_ratio_3_2,3), round(err_ratio_3_2/ratio_3_2*100.,2), '%'
    if not newscales: # FIXME #
        print 'PDFs up =', round(err_pdf_3_2_up,3), round(err_pdf_3_2_up/ratio_3_2*100.,2), '%'
        print 'PDFs down =', round(err_pdf_3_2_down,3), round(err_pdf_3_2_down/ratio_3_2*100.,2), '%'
    print 'extr up =', round(err_extr_3_2_up,3), round(err_extr_3_2_up/ratio_3_2*100.,2), '%'
    print 'extr down =', round(err_extr_3_2_down,3), round(err_extr_3_2_down/ratio_3_2*100.,2), '%'
    if do_scale_variations:
        print 'scale up =', round(err_scale_3_2_up,3), round(err_scale_3_2_up/ratio_3_2*100.,2), '%'
        print 'scale down =', round(err_scale_3_2_down,3), round(err_scale_3_2_down/ratio_3_2*100.,2), '%'
    print 'total =', round(.5*(err_3_2_up+err_3_2_down)/ratio_3_2*100.,2), '%'
    print '\n'

    print '\n'
    print 'uncertainties ratio_4_2:\n'
    print 'experimental =', round(err_ratio_4_2,3), round(err_ratio_4_2/ratio_4_2*100.,2), '%'
    if not newscales: # FIXME #
        print 'PDFs up =', round(err_pdf_4_2_up,3), round(err_pdf_4_2_up/ratio_4_2*100.,2), '%'
        print 'PDFs down =', round(err_pdf_4_2_down,3), round(err_pdf_4_2_down/ratio_4_2*100.,2), '%'
    print 'extr up =', round(err_extr_4_2_up,3), round(err_extr_4_2_up/ratio_4_2*100.,2), '%'
    print 'extr down =', round(err_extr_4_2_down,3), round(err_extr_4_2_down/ratio_4_2*100.,2), '%'
    if do_scale_variations:
        print 'scale up =', round(err_scale_4_2_up,3), round(err_scale_4_2_up/ratio_4_2*100.,2), '%'
        print 'scale down =', round(err_scale_4_2_down,3), round(err_scale_4_2_down/ratio_4_2*100.,2), '%'
    print 'total =', round(.5*(err_4_2_up+err_4_2_down)/ratio_4_2*100.,2), '%'
    print '\n'

    print 'results:\n'
    print 'ratio_1_2 =', round(ratio_1_2,3), '+' , round(err_1_2_up,3), '-' , round(err_1_2_down,3)
    print 'ratio_3_2 =', round(ratio_3_2,3), '+' , round(err_3_2_up,3), '-' , round(err_3_2_down,3)
    print 'ratio_4_2 =', round(ratio_4_2,3), '+' , round(err_4_2_up,3), '-' , round(err_4_2_down,3) 
    print
    
    return [err_1_2_up, err_1_2_down, err_3_2_up, err_3_2_down, err_4_2_up, err_4_2_down]

################################

# main program

################################


def execute():

    global doingToys
    doingToys = False

    setMasses()
    if ntoys>0 and replace_corr: estimateMassCorrelations()

    mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0 , 0)
    mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0 , 0)
    mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0 , 0)
    mass_and_err_4 = getMassAndError(4, 'nominal', 'nominal', 0 , 0 , 0)

    ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0], mass_and_err_4[0],
                                mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1], mass_and_err_4[1])

    ratio_1_2 = ratios_and_errs[0]
    ratio_3_2 = ratios_and_errs[1]
    ratio_4_2 = ratios_and_errs[2]

    # FIXME #
    if not newscales: pdf_errors = getPDFUncertainties(ratio_1_2, ratio_3_2, ratio_4_2)
    else: pdf_errors = []
    
    extr_errors = getExtrapolationUncertainties(ratio_1_2, ratio_3_2, ratio_4_2)

    if do_scale_variations: scale_errors = getScaleUncertainties(ratio_1_2, ratio_3_2, ratio_4_2)
    else : scale_errors = []

    if not newscales: # FIXME #
        makeMassPlots(mass_and_err_1[0], mass_and_err_1[1], mass_and_err_1[2],
                      mass_and_err_2[0], mass_and_err_2[1], mass_and_err_2[2],
                      mass_and_err_3[0], mass_and_err_3[1], mass_and_err_3[2],
                      mass_and_err_4[0], mass_and_err_4[1], mass_and_err_4[2],
        )

    err_ratio_1_2 = ratios_and_errs[3]
    err_ratio_3_2 = ratios_and_errs[4]
    err_ratio_4_2 = ratios_and_errs[5]

    tot_err = getTotalError (ratios_and_errs, pdf_errors, extr_errors, scale_errors)
    err_1_2_up = tot_err[0]
    err_1_2_down = tot_err[1]
    err_3_2_up = tot_err[2]
    err_3_2_down = tot_err[3]
    err_4_2_up = tot_err[4]
    err_4_2_down = tot_err[5]
    
    makeRatioPlots (mass_and_err_2[0], ratio_1_2, ratio_3_2, ratio_4_2, err_1_2_up, err_1_2_down, err_3_2_up, err_3_2_down, err_4_2_up, err_4_2_down, mass_and_err_2[2])
    # makeChi2Test (mass_and_err_2[0], ratio_1_2, ratio_3_2, ratio_4_2, err_1_2_up, err_1_2_down, err_3_2_up, err_3_2_down, err_4_2_up, err_4_2_down)
    if estimate_contribs : estimateSystContributions (ratio_1_2,ratio_3_2,ratio_4_2)
    if estimate_significance : makeChi2Significance (mass_and_err_2[0], ratio_1_2, ratio_3_2, ratio_4_2, err_ratio_1_2, err_ratio_3_2, err_ratio_4_2)

    return


################################

# execute main program

################################


execute()









