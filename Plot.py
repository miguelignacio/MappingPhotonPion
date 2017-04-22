import ROOT
from ROOT import TFile

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib import colors as mcolors

from AtlasCommonUtils import *
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['axes.labelsize']  = 15
#mpl.rcParams['legend.fontsize'] = 0.5
#mpl.rcParams['legend.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3.8
mpl.rcParams['legend.numpoints'] =1
mpl.rcParams['legend.fontsize'] = 'large'

SetAtlasStyle()
def SharkFin(x, a, b):
    #return 2*(a-b)*pionPT
    if(x <a): return 0
    elif(x >= a and x <=b): return (1-a/x)
    elif(x >b): return (b-a)/x
    return 0

def FromTH1toArray(histo):
  print ' Transforming histogram '
  # Calculate how many y bins we have
  xaxis = histo.GetXaxis()
  nx = xaxis.GetNbins()
  xbins = range(1, nx+1) # Removes under/overflow
  d = {}
  d["x_center"] = np.array([xaxis.GetBinCenter(i) for i in xbins])
  d["x_low"]    = np.array([xaxis.GetBinLowEdge(i) for i in xbins])
  d["x_up"]     = np.array([xaxis.GetBinUpEdge(i) for i in xbins])
  d["y"]        = np.array([histo.GetBinContent(i) for i in xbins])
  d['dy']       = np.array([histo.GetBinError(i) for i in xbins])
  return d

fin = ROOT.TFile("fout.root","READ")
h = fin.Get("Hsparse")


c = ROOT.TCanvas()

h.Projection(1,2).Draw("colz")
c.SaveAs("correlation.png")
c.Clear()

h.GetAxis(2).SetRange(11,14)
mass = h.Projection(3)
mass.Fit("gaus")
fit = mass.GetFunction("gaus")
fit.SetLineColor(2)
mass.Draw()
mass.GetXaxis().SetNdivisions(6)
mass.SetTitle(' ; $M_{\gamma\gamma}$ [GeV]')
fit.Draw("same")
ROOT.gPad.SetLogy(1)
c.SaveAs("masa.png")
h.GetAxis(2).SetRange(0,100)

den = h.Projection(0)



h.GetAxis(1).SetRange(11,14)
num = h.Projection(0)
h_ratio = num
h_ratio.Divide(den)

ratio = FromTH1toArray(h_ratio)

h.GetAxis(1).SetRange(0,100)
h.GetAxis(2).SetRange(11,14)
num = h.Projection(0)
h_ratio = num
h_ratio.Divide(den)

ratio2 = FromTH1toArray(h_ratio)


#f, axes = plt.subplots(1,4, sharex=True, sharey=True)
xfit = np.linspace(0.0, 20.0, num=1000, endpoint=True)
model = [SharkFin(x,5,7) for x in xfit]
x  = ratio["x_center"]
y  = ratio["y"]
dy = ratio["dy"]
plt.errorbar(x,y,yerr=dy, fmt='o', label='MC', color='crimson')                                
x  = ratio2["x_center"]
y  = ratio2["y"]
dy = ratio2["dy"]

#plt.errorbar(x,y,yerr=dy, fmt='--', label='MC, Smeared', color='blue')
plt.plot(xfit, model, linestyle='--', label ='Analytic', color = 'orange')

axes = plt.gca()
axes.set_xlim([4.0,17])
#axes.set_ylim([ymin,ymax])
plt.legend(loc='best')
plt.xlabel(r' $\pi^{0} p_{\mathrm{T}}$ [GeV]')
plt.ylabel('P')
plt.tight_layout()
#plt.show()
plt.savefig('Sharkfin_3.png')

plt.clf()
plt.errorbar(x,y,yerr=dy, fmt='--', label='MC, Smeared', color='blue')
#plt.show()
plt.savefig('Sharkfin_2.png')

print 'ola'

  


