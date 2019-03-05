import ROOT
from ROOT import *
import os
import rundec


# input values (global variables)

estimate_contribs = True
ntoys = 0

replace_corr = False

n_mttbins = 4

mass_min = 152.0
mass_max = 172.0
mass_fine_min = 152.0
mass_fine_max = 172.0

mu_1 = 383.95
mu_2 = 476.19
mu_3 = 644.27

MZ = 91.1876
asMZ = 0.11905
err_as = 0.0011
nflav = 5
nloops = 2


xsec_1 = 248.98
err_xsec_1_up = 5.759
err_xsec_1_down = 5.333
extr_1_up = [ -0.14 , 2.545 , -0.118 , -0.049 , 2.107 , -1.832 , 1.504 , ]
extr_1_down = [ 0.663 , -0.274 , 0.448 , 0.262 , -1.422 , 0.519 , -3.826 , ]
err_xsec_toys_1 = 3.818

extr_name = [ 'ISR_scale' , 'FSR_scale' , 'ME_scale' , 'UE_tune' , 'PDF' , 'top_pt' , 'mt' , ]

syst_names = [ 'b-tag_b_fragm' , 'b-tag_b_templ_corr' , 'b-tag_JP_corr' , 'b-tag_c_fragm' , 'b-tag_D_tp_mu_X_BR' ,
               'b-tag_gluon_split' , 'b-tag_JES' , 'b-tag_away_jet_tag' , 'b-tag_Ks0_Lambda_prod' , 'b-tag_l/c_ratio' ,
               'b-tag_LT_others' , 'b-tag_muon_DR' , 'b-tag_muon_pT' , 'b-tag_muon_pT_rel' , 'b-tag_sample_dep' ,
               'b-tag_stat' , 'b-mistag' , 'DY_ME_scale' , 'Electron_ER_phi' , 'Electron_ER_rho' ,
               'Electron_energy_scale' , 'Electron_ID_isolation' , 'Jet_energy_resolution' , 'JES_MPF' , 'JES_Absolute_Scale' ,
               'JES_Absolute_Stat' , 'JES_Fragmentation' , 'JES_Pileup_Data/MC' , 'JES_Pileup_pT_BB' , 'JES_PileUp_pT_EC1' ,
               'JES_PileUp_pT_Ref' , 'JES_Relative_Balance' , 'JES_Intercalibration' , 'JES_Relative_JER_EC1' , 'JES_Relative_pT_BB' ,
               'JES_Relative_pT_EC1' , 'JES_Relative_Stat_EC' , 'JES_Relative_Stat_FSR' , 'JES_Single_pion_ECAL' , 'JES_Single_pion_HCAL' ,
               'JES_Time_pT_eta' , 'Muon_energy_scale' , 'Muon_ID_isolation' , 'Pile-up' , 'tW_FSR_scale' ,
               'tW_ISR_scale' , 'tW_ME_scale' , 'top_mass' , 'Top_pT' , 'Trigger' ,
               'B-hadron_BR' , 'CR_ERD_on' , 'CR_Gluon_move' , 'CR_QCD-inspired' , 'fragm_Peterson' ,
               'fragmentation' , 'ttbar_FSR_scale' , 'NLO_generator' , 'ttbar_ISR_scale' , 'ME/PS_matching' ,
               'ttbar_ME_scale' , 'UE_tune' , 'PDF_10' , 'PDF_11' , 'PDF_12' ,
               'PDF_13' , 'PDF_14' , 'PDF_15' , 'PDF_16' , 'PDF_17' ,
               'PDF_18' , 'PDF_19' , 'PDF_1' , 'PDF_20' , 'PDF_21' ,
               'PDF_22' , 'PDF_23' , 'PDF_24' , 'PDF_25' , 'PDF_26' ,
               'PDF_27' , 'PDF_28' , 'PDF_2' , 'PDF_3' , 'PDF_4' ,
               'PDF_5' , 'PDF_6' , 'PDF_7' , 'PDF_8' , 'PDF_9' ,
               'JES_Flavor' , 'tW_background' , 'W+jets_background' , 'Diboson_background' , 'ttbar_background' ,
               'DY_bg_0_b-jets' , 'DY_bg_1_b-jets' , 'DY_bg_2_b-jets' , 'Luminosity' , ]

contribs_1 = [ -0.015 , -0.062 , -0.03 , 0.0 , -0.054 , -0.05 , -0.019 , -0.02 , -0.001 , -0.026 ,
               -0.184 , 0.006 , -0.014 , -0.043 , -0.236 , -0.133 , 0.0 , -0.032 , 0.275 , -0.071 ,
               -0.009 , -1.853 , 0.307 , 0.033 , -0.023 , -0.101 , -0.211 , -0.159 , -0.5 , -0.04 ,
               0.028 , 0.135 , -0.02 , -0.7 , -0.373 , -0.109 , 0.152 , 0.048 , -0.16 , -0.388 ,
               -0.257 , -0.012 , -1.357 , 0.595 , 0.032 , -0.063 , -0.113 , 0.314 , -3.019 , -0.357 ,
               0.102 , -0.859 , 0.683 , 0.022 , 0.662 , -0.603 , 0.311 , 0.0 , 0.172 , 0.189 ,
               -0.605 , 0.159 , 0.749 , 0.225 , -0.404 , -0.347 , -0.302 , -0.213 , 0.077 , -0.099 ,
               0.645 , 0.759 , -0.03 , 0.25 , 0.016 , 0.15 , -0.179 , -0.206 , 0.353 , 0.006 ,
               0.245 , 0.009 , 0.518 , 0.302 , 0.246 , -0.873 , 0.151 , -1.155 , -0.297 , 0.208 ,
               0.042 , -0.943 , -0.065 , 0.093 , -0.193 , 0.18 , -0.122 , -0.001 , -2.696 , ]

xsec_2 = 310.3
err_xsec_2_up = 4.779
err_xsec_2_down = 4.469
extr_2_up = [ -0.131 , 1.717 , -0.346 , -0.111 , 1.33 , 0.021 , -2.224 , ]
extr_2_down = [ 0.459 , -0.166 , 0.707 , 0.029 , -0.906 , -0.071 , 0.809 , ]
err_xsec_toys_2 = 4.918

contribs_2 = [ 0.148 , 0.052 , 0.012 , 0.0 , 0.029 , -0.083 , -0.028 , -0.151 , 0.002 , 0.028 ,
               0.174 , 0.004 , 0.024 , -0.026 , 0.081 , 0.002 , 0.0 , -0.127 , -0.241 , 0.032 ,
               0.085 , -1.874 , -0.448 , -0.226 , -0.237 , -0.141 , 0.06 , -0.128 , -0.253 , -0.186 ,
               -0.301 , -0.367 , 0.146 , 0.694 , 0.087 , 0.133 , -0.04 , 0.155 , -0.297 , 0.311 ,
               -0.013 , -0.203 , -1.31 , 0.271 , 0.091 , 0.05 , -0.002 , -0.062 , 0.104 , -0.358 ,
               0.134 , 0.38 , -0.335 , 0.595 , 0.504 , -0.745 , 1.488 , 0.0 , -0.655 , -0.849 ,
               -1.37 , -0.17 , 0.579 , 0.199 , -0.324 , 0.057 , -0.168 , -0.126 , 0.084 , -0.206 ,
               0.214 , 0.494 , -0.042 , 0.201 , -0.124 , -0.117 , -0.127 , 0.015 , 0.126 , 0.014 ,
               0.13 , -0.032 , 0.228 , 0.193 , 0.186 , -0.541 , 0.216 , -0.724 , -0.114 , -0.002 ,
               -0.035 , -0.875 , -0.059 , -0.074 , 0.109 , 0.771 , 0.027 , -0.005 , -2.644 , ]

xsec_3 = 186.91
err_xsec_3_up = 5.527
err_xsec_3_down = 5.176
extr_3_up = [ -0.187 , 1.609 , 0.652 , 0.147 , 1.031 , 0.689 , -3.36 , ]
extr_3_down = [ 0.019 , -0.07 , -0.125 , -0.091 , -0.869 , -2.306 , 1.22 , ]
err_xsec_toys_3 = 2.68

contribs_3 = [ 0.123 , 0.081 , 0.018 , 0.001 , 0.077 , -0.175 , -0.107 , -0.229 , 0.002 , 0.052 ,
               0.249 , -0.013 , 0.029 , 0.001 , 0.168 , 0.061 , 0.0 , -0.123 , 0.598 , -0.138 ,
               -0.318 , -1.826 , 0.091 , 0.206 , 0.295 , -0.364 , -0.76 , 0.258 , -0.364 , -0.38 ,
               -0.177 , 0.255 , -0.3 , -0.563 , 0.037 , -0.389 , -0.076 , 0.075 , -0.29 , -0.068 ,
               -0.008 , 0.113 , -1.172 , 0.181 , 0.104 , -0.022 , -0.143 , 0.371 , 3.186 , -0.375 ,
               -0.029 , 0.411 , -0.1 , 0.314 , 0.349 , -1.055 , 0.626 , 0.0 , 0.049 , 0.194 ,
               -0.079 , -0.084 , 0.266 , 0.146 , -0.166 , 0.591 , 0.033 , 0.026 , 0.093 , -0.322 ,
               -0.398 , 0.047 , -0.049 , 0.054 , -0.348 , -0.552 , 0.018 , 0.338 , -0.229 , 0.028 ,
               -0.039 , -0.057 , -0.286 , 0.063 , 0.007 , -0.04 , 0.188 , -0.054 , 0.171 , -0.246 ,
               -0.129 , -0.523 , -0.065 , 0.014 , 0.139 , 0.619 , 0.09 , 0.013 , -2.594 , ]

corr_1_2 = 0.57
corr_3_2 = 0.47

# end global variables


################################

# "setMasses" initializes which masses should be considered

################################

def setMasses():

    global mass_v
    mass_v = []
    i_mass = mass_max
    while i_mass >= mass_min :
        mass_v.append(i_mass)
        if i_mass <= mass_fine_max and i_mass > mass_fine_min+.1 : i_mass-= 0.2
        else : i_mass-= 1.0
    return


################################

# "mtmt2mtmu" converts mt(mt) to mt(mu) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmt2mtmu(mt, mu):

    crd = rundec.CRunDec()

    asmt = crd.AlphasExact(asMZ, MZ, mt, nflav, nloops)
    asmu = crd.AlphasExact(asMZ, MZ, mu, nflav, nloops)
    
    mtmu = crd.mMS2mMS(mt, asmt, asmu, nflav, nloops)
    
    return mtmu



################################

# "mtmu2mtmu" converts mt(mu1) to mt(mu2) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmu2mtmu (mtmu1, mu1, mu2, var_as):

    alphaS = 0
    if var_as == 'nominal' : alphaS = asMZ
    elif var_as == 'up' : alphaS = asMZ+err_as
    elif var_as == 'down' : alphaS = asMZ-err_as
    else :
        print 'ERROR!'
        return 0

    crd = rundec.CRunDec()
    
    asmu1 = crd.AlphasExact(alphaS, MZ, mu1, nflav, nloops)
    asmu2 = crd.AlphasExact(alphaS, MZ, mu2, nflav, nloops)
    
    mtmu2 = crd.mMS2mMS(mtmu1, asmu1, asmu2, nflav, nloops)
    
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
    elif n_mttbins == 3:
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

    graph = ROOT.TGraph()
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
            xsec_err = ((xsec_1/100.*(err_xsec_1_up+err_xsec_1_down)/2.)**2 + err_xsec_toys_1**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_1_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_1_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_1_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_1_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + contribs_1[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - contribs_1[abs(contrib)-1]/100.
                
        elif mttbin == 2 : 
            xsec_exp = xsec_2 
            xsec_err = ((xsec_2/100.*(err_xsec_2_up+err_xsec_2_down)/2.)**2 + err_xsec_toys_2**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_2_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_2_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_2_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_2_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + contribs_2[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - contribs_2[abs(contrib)-1]/100.

        else : # mttbin = 3
            xsec_exp = xsec_3 
            xsec_err = ((xsec_3/100.*(err_xsec_3_up+err_xsec_3_down)/2.)**2 + err_xsec_toys_3**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_3_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_3_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_3_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_3_down[abs(extrapol)-1]/100.
            if contrib !=0 :
                if contrib > 0 : xsec_exp *= 1 + contribs_3[abs(contrib)-1]/100.
                else : xsec_exp *= 1 - contribs_3[abs(contrib)-1]/100.


        # chi2 = abs(xsec_th-xsec_exp)/xsec_err
        chi2 = (xsec_th-xsec_exp)/xsec_err
        fact_A = 0
        if mttbin == 1 : fact_A = (err_xsec_1_up-err_xsec_1_down)/(err_xsec_1_up+err_xsec_1_down)
        elif mttbin == 2 : fact_A = (err_xsec_2_up-err_xsec_2_down)/(err_xsec_2_up+err_xsec_2_down)
        else : fact_A = (err_xsec_3_up-err_xsec_3_down)/(err_xsec_3_up+err_xsec_3_down)
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

        mass_v_sel.append(mass)
        chi2_v.append(chi2)
            
    for i in range(0,len(chi2_v)):
        graph.SetPoint(i,mass_v_sel[i],chi2_v[i])

    funct = TF1('funct','pol4(0)',mass_min,mass_max)
    funct.SetParameter(0,1080)
    funct.SetParameter(1,-12)
    funct.SetParameter(2,0.033)
    # graph.Fit(funct,'q','',163,166)
    graph.Fit(funct,'q','',mass_min+0.1,mass_max-0.1)

    line_up = ROOT.TLine(mass_min,1.,mass_max,1.)
    line_down = ROOT.TLine(mass_min,-1.,mass_max,-1.)
    line_central = ROOT.TLine(mass_min,0.,mass_max,0.)
    line_up.SetLineColor(kGreen+2)
    line_down.SetLineColor(kGreen+2)
    line_central.SetLineColor(kRed)

    outdir = 'plots_chi2'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c = ROOT.TCanvas()
    c.SetGrid()
    graph.Draw('ap')
    graph.SetMarkerStyle(7)
    graph.SetTitle('signed #chi^{2} vs top mass in MCFM prediction; m_{t}(m_{t}) [GeV]; data-theory #chi^{2}')
    line_up.Draw('same')
    line_down.Draw('same')
    line_central.Draw('same')
    if not doingToys:
        c.Print(outdir+'/test_mtt'+str(mttbin)+'_mur_'+murscale+'_muf_'+mufscale+'_extrapol'+str(extrapol)+'_pdf'+str(pdfmember)+'_contrib'+str(contrib)+'.png')

    fitted_mass = funct.GetX(0)
    fitted_mass_up = funct.GetX(-1)
    fitted_mass_down = funct.GetX(1)

    #now evolve masses
    mu = 0
    if mttbin == 1 : mu = mu_1
    elif mttbin == 2 : mu = mu_2
    else : mu = mu_3 # mttbin = 3

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0 and contrib==0 and not doingToys:
        print
        print 'mt(mt) bin', mttbin,'=', round(fitted_mass,2), round(fitted_mass_up-fitted_mass,2), round(fitted_mass-fitted_mass_down,2)
    
    fitted_mass = mtmt2mtmu(fitted_mass, mu)
    fitted_mass_up = mtmt2mtmu(fitted_mass_up, mu)
    fitted_mass_down = mtmt2mtmu(fitted_mass_down, mu)

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0 and contrib==0 and not doingToys:
        print 'mt(mu) bin', mttbin,'=', round(fitted_mass,2), round(fitted_mass_up-fitted_mass,2), round(fitted_mass-fitted_mass_down,2)
        
    fitted_mass_err = (fitted_mass_up - fitted_mass_down)/2

    return [fitted_mass, fitted_mass_err]


################################

# "getRatios" calculates all the ratios and their experimental errors
# taking the correlations into account

################################

def getRatios(mass_1, mass_2, mass_3, err_1, err_2, err_3):

    ratio_1_2 = mass_1/mass_2
    err_ratio_1_2 = (err_1/mass_1)**2 + (err_2/mass_2)**2 - 2*corr_1_2*(err_1/mass_1)*(err_2/mass_2)
    err_ratio_1_2 = err_ratio_1_2**.5
    err_ratio_1_2*= ratio_1_2

    ratio_3_2 = mass_3/mass_2
    err_ratio_3_2 = (err_3/mass_3)**2 + (err_2/mass_2)**2 - 2*corr_3_2*(err_3/mass_3)*(err_2/mass_2)
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
    for extr in range(-len(extr_1_up),len(extr_1_up)+1) :
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

        name = extr_name[abs(extr)-1]
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
        ratio = mtmu2mtmu(mass_2, mu_2, scale, 'nominal')/mass_2
        out.write(str(scale)+'\t'+str(ratio)+'\n')
        ratio_up = mtmu2mtmu(mass_2, mu_2, scale, 'up')/mass_2
        ratio_down = mtmu2mtmu(mass_2, mu_2, scale, 'down')/mass_2
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

# "makePlots" produces the plots for the running of mt

################################


def makePlots (mass_2, ratio_12, ratio_32, err_12, err_32):

    graph = TGraphErrors(3)

    graph.SetPoint(0,mu_1,ratio_12)
    graph.SetPointError(0,0,err_12)

    graph.SetPoint(1,mu_2,1)
    graph.SetPointError(1,0,0)

    graph.SetPoint(2,mu_3,ratio_32)
    graph.SetPointError(2,0,err_32)

    theoryFileName = 'theory_prediction'
    l = makeTheoryPrediction(theoryFileName,mass_2)
    th = TGraph(theoryFileName+'.txt')
    th.SetLineColor(kRed);
    th.SetLineWidth(2);

    scales = l[0]
    r_up = l[1]
    r_down = l[2]
    
    th_band = TGraph (2*th.GetN())

    for i in range(0, th.GetN()):
        th_band.SetPoint(i,scales[i],max(r_up[i],r_down[i]))
        th_band.SetPoint(th.GetN()+i,scales[th.GetN()-i-1],min(r_up[th.GetN()-i-1],r_down[th.GetN()-i-1]));

    th_band.SetFillStyle(3001);
    th_band.SetFillColor(kRed);
    
    leg = TLegend(.15,.2,.6,.32)
    leg.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg.AddEntry(th,'RunDec @ 2 loops (5 flav.)','l')

    latexLabel1 = TLatex();
    latexLabel1.SetTextSize(0.04);
    latexLabel1.SetTextFont(42)
    latexLabel1.SetNDC();
    
    graph.GetXaxis().SetTitle('median m_{t#bar{t}} [GeV]')
    graph.GetYaxis().SetTitle('m_{t}(m_{t#bar{t}}) / m_{t}('+str(int(mu_2))+' GeV)')
    graph.SetTitle('running of m_{t}(#mu) as a function of #mu=m_{t#bar{t}}')

    c = TCanvas()
    graph.SetMarkerStyle(8)
    graph.Draw('ap')
    th.Draw("L same");
    leg.Draw('same')
    th_band.Draw('f same')
    latexLabel1.DrawLatex(0.59, 0.78, "ABMP16_5_nlo PDF set");

    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/running.png')
    c.SaveAs(outdir+'/running.pdf')
    c.SaveAs(outdir+'/running.root')

    mtmt = 164.02
    mtmt_err = 1.58

    mtmu = mtmt2mtmu (mtmt, mu_2)

    gr_add = TGraphErrors()

    gr_add.SetPoint(0,mtmt,mtmt/mtmu)
    gr_add.SetPointError(0,0,mtmt_err/mtmu)

    l = makeAdditionalTheoryPrediction(mtmt, mtmt_err, mtmu)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]

    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1]);

    gr_band.SetFillStyle(3002);
    gr_band.SetFillColor(kRed+1);

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(#mu) / m_{t}('+str(int(mu_2))+' GeV)')
    gr_band.SetTitle('running of m_{t}(#mu) as a function of #mu')

    gr_add.SetMarkerStyle(8)
    gr_add.SetMarkerColor(kBlue)
    gr_add.SetLineColor(kBlue)

    gr_band.GetYaxis().SetRangeUser(0.9,1.1)

    leg2 = leg.Clone()
    leg2.Clear()
    leg2.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg2.AddEntry(gr_add,'Hathor @NLO from incl. #sigma_{t#bar{t}} (same data)')
    leg2.AddEntry(gr_band,'RunDec @ 2 loops (5 flav.)','f')

    g = graph.Clone()
    g.RemovePoint(1)
    g1 = TGraph()
    g1.SetPoint(0,mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c.Clear()
    gr_band.Draw('af')
    gr_add.Draw('p same')
    g.Draw('p same')
    g1.Draw('p same')
    # th.Draw("L same");
    leg2.Draw('same')
    # th_band.Draw('f same')

    c.SaveAs(outdir+'/test_incl.png')
    c.SaveAs(outdir+'/test_incl.pdf')
    c.SaveAs(outdir+'/test_incl.root')

    return



################################

# "makeChi2Test" calculates the significance of the obseved running
#  using a chi2 test

################################


def makeChi2Test(mass2, ratio12, ratio32, err12, err32):

    th_ratio12 = mtmu2mtmu(mass2, mu_2, mu_1, 'nominal')/mass2
    th_ratio32 = mtmu2mtmu(mass2, mu_2, mu_3, 'nominal')/mass2

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

    for contrib in range(1,len(contribs_1)+1) :

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
            print syst_names[contrib-1], round(contr_1_2,2), round(contr_3_2,2)


    print
    
    return



################################

# "throwToyCrossSections" throws toy cross sections taking the correlations into account
#  the results are used to estimate the correlation between the fitted masses

################################

def throwToyCrossSections(r):

    global xsec_1, xsec_2, xsec_3

    orig_xsec_1 = xsec_1
    orig_xsec_2 = xsec_2
    orig_xsec_3 = xsec_3

    err_1 = xsec_1*(err_xsec_1_up+err_xsec_1_down)/2./100.
    err_2 = xsec_2*(err_xsec_2_up+err_xsec_2_down)/2./100.
    err_3 = xsec_3*(err_xsec_3_up+err_xsec_3_down)/2./100.

    xsec_2 = r.Gaus(xsec_2,err_2)
    xsec_1 = r.Gaus(xsec_1,err_1*(1-corr_1_2*corr_1_2)**.5) + corr_1_2*err_1/err_2*(xsec_2-orig_xsec_2)
    xsec_3 = r.Gaus(xsec_3,err_3*(1-corr_3_2*corr_3_2)**.5) + corr_3_2*err_3/err_2*(xsec_2-orig_xsec_2)
    
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

    r_corr = TH2D('r_corr','r_corr',100,0.9,1.15,100,0.75,1.15)
    
    
    for i in range(1,ntoys+1):
        if i%1000==0 or i==1: print 'toy n.', i
        throwToyCrossSections(r)
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
        global corr_1_2, corr_3_2
        corr_1_2 = m12.GetCorrelationFactor()
        corr_3_2 = m32.GetCorrelationFactor()
    
    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c = TCanvas()
    m12.SetTitle('correlation = '+str(round(m12.GetCorrelationFactor(),2)))
    m12.GetXaxis().SetTitle('mt ('+str(int(mu_1))+' GeV) [GeV]')
    m12.GetYaxis().SetTitle('mt ('+str(int(mu_2))+' GeV) [GeV]')
    m12.DrawNormalized('colz')
    c.SaveAs(outdir+'/mass_corr_1_2.png')
    c.SaveAs(outdir+'/mass_corr_1_2.pdf')
    c.SaveAs(outdir+'/mass_corr_1_2.root')
    c.Clear()
    m32.SetTitle('correlation = '+str(round(m32.GetCorrelationFactor(),2)))
    m32.GetXaxis().SetTitle('mt ('+str(int(mu_3))+' GeV) [GeV]')
    m32.GetYaxis().SetTitle('mt ('+str(int(mu_2))+' GeV) [GeV]')
    m32.DrawNormalized('colz')
    c.SaveAs(outdir+'/mass_corr_3_2.png')
    c.SaveAs(outdir+'/mass_corr_3_2.pdf')
    c.SaveAs(outdir+'/mass_corr_3_2.root')
    c.Clear()
    r_corr.SetTitle('correlation = '+str(round(r_corr.GetCorrelationFactor(),2)))
    r_corr.GetXaxis().SetTitle('mt ('+str(int(mu_1))+' GeV) / mt ('+str(int(mu_2))+' GeV)')
    r_corr.GetYaxis().SetTitle('mt ('+str(int(mu_3))+' GeV) / mt ('+str(int(mu_2))+' GeV)')
    r_corr.DrawNormalized('colz')
    c.SaveAs(outdir+'/ratio_corr.png')
    c.SaveAs(outdir+'/ratio_corr.pdf')
    c.SaveAs(outdir+'/ratio_corr.root')

    print
    if replace_corr:
        print 'new corr_1_2 =', round(m12.GetCorrelationFactor(),3)
        print 'new corr_3_2 =', round(m32.GetCorrelationFactor(),3)
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

    makePlots (mass_and_err_2[0], ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2)
    makeChi2Test (mass_and_err_2[0], ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2)
    if estimate_contribs : estimateSystContributions (ratio_1_2,ratio_3_2)
    
    return


################################

# execute main program

################################


execute()









