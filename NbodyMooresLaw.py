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

# N-body simulation data
MIN = 60
HR = 60 * MIN
DAY = 24 * HR
MONTH = 30 * DAY
papers = {
    1950:((40/365)/(2*MIN),0,"Ekert '52",0,"EBC52"),
    1965:((1500)/(HR),0,"Cohen & Hubbard '65",0,"CH65"),
    1978:((40/365)/(10),0,"TRS-80",0,"TRS-80"),
    1984.1:(5e6/(4*HR),0,"Kinoshita & Nakai'84",1,"KN84"),
    1986:(60*1e8/(365*DAY),0,"Applegate+ '86",1,"A+86"),
    1991.1:(3e6/(2*MONTH),1,"Quinn+ '91",0,"QDT91"),
    1991:(1e9/(14*DAY),0,"Wisdom & Holman '91",0,"WH91"),
    2008:(20e9 / (6*MONTH),1,"Batygin & Laughlin '08",0,"BL08"),
    2009:(5e9 / (2500*HR),1,"Laskar & Gastineau '09",0,"LG09"),
    2020:(5e9*96 / (6*12*MONTH),1,"Brown & Rein '20",0,"BR20"),
    2023:(2 * 2750 * 5e9 / (2.5e6*HR) ,1,"Abbot+ '23",0,"A+23"),
    2023.1:(1e9 / (DAY) ,1,"Javaheri+ '23",0,"JRT23"),
}

# CPU clock speed data
with open("./frequency.dat","r") as fi:
    lines = fi.readlines()
freq_data=np.array([list(map(float,l.split())) for l in lines])
yr,freq=freq_data.T

# Time re-scaling for simulations of outer SS only
#   -- Note- should probably include factor of something like (8/4)**2 to be more fair
OUTER_TO_INNER_RESCALE = (5.2/0.387)**(1.5)


# the plot
fig=plt.figure(figsize=(10,7))
ax = plt.gca()
plt.tick_params(labelsize=16,size=8,direction='in')
plt.tick_params(size=6,which='minor',direction='in')
TO_MYR_PER_MONTH = MONTH/1e6
for year,data in papers.items():
    rate,inner,label,hardware,shortname=data
    if inner:
        x,y=year,TO_MYR_PER_MONTH * rate
        if year==2023:
            plt.scatter(x,y,color='k',zorder=99,marker='*',s=150,label='Simulation Efficiency [Myr/CPU month]')
        else:
            plt.scatter(x,y,color='k',zorder=99,marker='*',s=150)
        plt.text(x+0.2,y*1.25,shortname,ha='center',fontsize=12)
    else:
        x,y = year,TO_MYR_PER_MONTH * rate / OUTER_TO_INNER_RESCALE 
        plt.scatter(x,y,color='k',marker='*',s=150)
        plt.text(x+0.2,y*1.25,shortname,ha='center',zorder=99,fontsize=12)
# known planets
ax.plot(
    np.concatenate(([1940],np.sort(disc_years))),
    np.arange(8,1+8+len(disc_years))
    ,'k-',
    lw=3,
    label='Known Planets'
)
# CPU clock speed
ax.plot(yr,freq,'b.',ms=3,zorder=0,label='CPU clock rates [MHz]')

# Cosmetics
plt.yscale('log')
plt.ylabel("$N$",fontsize=16)
plt.xlim(1945,2035)
plt.xlabel("Year",fontsize=16)
plt.tick_params(labelsize=16,size=8,direction='in')
plt.tick_params(size=6,which='minor',direction='in')
plt.title("$N$-body simulations, computer efficiency, and planet discoveries",fontsize=16)
plt.legend(loc='lower right',fontsize=12)
plt.tight_layout()
outfile = "/Users/shadden/github_io_website/shadden.github.io/assets/images/Nbody-Moores-Law.png"
plt.savefig(outfile)
outfile="/Users/shadden/DropboxPersonal/Apps/Overleaf/UCSD_research_statement/Nbody-Moores-Law.pdf"
plt.savefig(outfile)
