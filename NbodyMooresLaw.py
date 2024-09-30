from matplotlib import pyplot as plt

import numpy as np
import pyvo as vo


# Exoplanet Archive Data
service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")
query="SELECT hostname,pl_name,disc_year FROM ps WHERE"
query+=" default_flag=1"
query
tab = service.search(query)
hosts = set(tab.getcolumn('hostname').data)
astro_table = tab.to_table()
disc_years=astro_table['disc_year'].data.data
MERC_PER_DIVIDED_BY_MOON_PER = 1/0.31
# N-body simulation data
MIN = 60
HR = 60 * MIN
DAY = 24 * HR
MONTH = 30 * DAY
S=8
WH512_eff = 1e9 / (8*HR) #5e9 / (DAY)
papers = {
    1952:((40/365)/(2*MIN),0,"Ekert '52",0,"EBC52",(S,S)),
    1965:((1500)/(HR),0,"Cohen & Hubbard '65",0,"CH65",(S,S)),
    1978:((40/365)/(10),0,"TRS-80",0,"TRS-80",(S,S)),
    1983:((MERC_PER_DIVIDED_BY_MOON_PER)*((10/9)**2)*(3002 - 1411)/(9*DAY),1,"Newhall, Standish, and Williams (1983)",0,"NSW83",(S,S)),
    1984.1:(5e6/(4*HR),0,"Kinoshita & Nakai'84",1,"KN84",(-5*S,S)),
    1986:(60*1e8/(365*DAY),0,"Applegate+ '86",1,"A+86",(S,S)),
    1991.1:(3e6/(2*MONTH),1,"Quinn+ '91",0,"QDT91",(S,S)),
    1991:(1e9/(14*DAY),0,"Wisdom & Holman '91",0,"WH91",(S,S)),
    2008:(20e9 / (6*MONTH),1,"Batygin & Laughlin '08",0,"BL08",(S,S)),
    2009:(5e9 / (2500*HR),1,"Laskar & Gastineau '09",0,"LG09",(S,S)),
    2020:(5e9*96 / (6*12*MONTH),1,"Brown & Rein '20",0,"BR20",(-5*S,S)),
    2023:(2 * 2750 * 5e9 / (2.5e6*HR) ,1,"Abbot+ '23",0,"A+23",(S,S)),
    2023.1:(WH512_eff,1,"Javaheri+ '23",0,"JRT23",(S,S)),
}

# CPU clock speed data
with open("./frequency.dat","r") as fi:
    lines = fi.readlines()
freq_data=np.array([list(map(float,l.split())) for l in lines])
yr,freq=freq_data.T

# Time re-scaling for simulations of outer SS only
#   -- Note- should probably include factor of something like (8/4)**2 to be more fair
OUTER_TO_INNER_RESCALE = (5.2/0.387)**(1.5) * (9/5)**2


# the plot
fig=plt.figure(figsize=(10,7))
ax = plt.gca()
plt.tick_params(labelsize=16,size=8,direction='in')
plt.tick_params(size=6,which='minor',direction='in')

# Cosmetics
plt.yscale('log')
ax.set_ylim(9e-6,10**(5.5))
ax.set_yticks(10.**np.arange(-5,6))
minor_tx = np.array([np.arange(2,10) * pow10 for pow10 in 10.**np.arange(-5,6)])
minor_tx = minor_tx.reshape(-1)
ax.set_yticks(minor_tx,minor=True)
plt.xlim(1945,2035)
plt.ylim(9e-06, 900000.0)
plt.tick_params(labelsize=20,size=8,direction='in')
plt.tick_params(size=6,which='minor',direction='in')
plt.ylabel("$N$",fontsize=24)
plt.xlabel("Year",fontsize=24)
plt.title("$N$-body simulations, CPU efficiency,\nand planet discoveries",fontsize=24)

plt.tight_layout()


# CPU clock speed
ax.plot(yr,freq,'bs',ms=8,zorder=0,label='CPU clock rates [MHz]')
plt.legend(loc='lower right',fontsize=16)
#plt.savefig("./Nbody_moores_2.pdf")

# N-body integrations
TO_MYR_PER_MONTH = MONTH/1e6
for year,data in papers.items():
    rate,inner,label,hardware,shortname,txtoffset=data
    if year==1952:
        lbl = "Simulation Efficiency\n[Myr/CPU month]"
    else:
        lbl=None
    if inner:
        x,y=year,TO_MYR_PER_MONTH * rate
        ax.scatter(x,y,marker='*',s=250,zorder=99,c='k',label=lbl)
        #plt.text(x,y,shortname,zorder=99,fontsize=12,backgroundcolor='black',color='white')
        ax.annotate(shortname,
                    xy=(x,y),
                    xytext=txtoffset,
                    textcoords='offset points',
                    color='black',
                    weight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=1)
                )
    else:
        x,y = year,TO_MYR_PER_MONTH * rate / OUTER_TO_INNER_RESCALE 
        #plt.text(x,y,shortname,ha='center',zorder=99,fontsize=12,backgroundcolor='black',color='white')
        ax.scatter(x,y,marker='*',s=250,zorder=99,c='k',label=lbl)
        ax.annotate(
            shortname,
            xy=(x,y),
            xytext=txtoffset, 
            textcoords='offset points',
            color='black',
            weight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=1)
            )
    print(label, np.log10(y))
plt.legend(loc='lower right',fontsize=16)
plt.savefig("./Nbody_moores_alt_1.pdf")

# known planets
ax.plot(
    np.concatenate(([1940],np.sort(disc_years))),
    np.arange(8,1+8+len(disc_years)),
    ls='-',
    color='orange',
    lw=6,
    label='Known Planets'
)
plt.legend()
plt.legend(loc='lower right',fontsize=16)
plt.savefig("./Nbody_moores_alt_2.pdf")



outfile = "/Users/shadden/github_io_website/shadden.github.io/assets/images/Nbody-Moores-Law.png"
plt.savefig(outfile)
# outfile="/Users/shadden/DropboxPersonal/Apps/Overleaf/UofI_research_statement_2023/Nbody-Moores-Law.pdf"
# plt.savefig(outfile)
# plt.show()
