"""
Convert BPASS spectra txt files to binary for input to FSPS
"""
import glob
import numpy as np
from hoki import load
from spectres import spectres


"""
User defined variables
"""
# Salpeter IMF
models = glob.glob('bpass_files/spectra-bin-imf135all*')
out_prep = 'bpass_v2.2_salpeter100'

# Chabrier IMF
# models = glob.glob('bpass_files/spectra-bin-imf_chab*')
# out_prep = 'bpass_v2.2.1_chabrier100'
"""
END
"""

output_temp = load.model_output(models[0])

ages = np.array([float(a) for a in output_temp.columns[1:]])
age_mask = (10**ages / 1e9) < 18 # Gyr
ages = ages[age_mask]

wl = output_temp['WL'].values

metallicities = np.array([None] * len(models))
spec = np.zeros((len(wl),len(ages),len(metallicities)))

for i,mod in enumerate(models):
    try:
        metallicities[i] = float('0.'+mod[-7:-4])
    except: # ...catch em# format
        metallicities[i] = 10**-float(mod[-5])


Z_idx = np.argsort(metallicities)
metallicities = metallicities[Z_idx].astype(float)

for i,(Z,mod) in enumerate(zip(metallicities,np.array(models)[Z_idx])):
    output = load.model_output(mod)
    for j,a in enumerate(ages):
        spec[:,j,i] = output[str(a)].values


# resample on original wavelength grid
wl_out = np.loadtxt('bpass.lambda')
wl_out = wl_out[1:-int(1e3)]

_spec = spectres(wl_out, wl, spec.reshape((len(spec),-1)).T)
_spec = _spec.T.reshape((len(wl_out),len(ages),len(metallicities)))

_spec *= 1e-6  # convert from 1e6 to 1 Msol
_spec *= wl_out[:,None,None]**2 * 1e-10 / 3e8 # convert from Lsol / AA -> Lsol / Hz

fname = out_prep+'.ssp.bin'
print("writing ",fname)
_spec.T.tofile(fname)

fname = out_prep+'.lambda'
print("writing ",fname)
np.savetxt(fname,wl_out)


"""
Save mass table
"""
files = glob.glob('bpass_files/starmass-bin-imf_chab100*')

z_arr = []
for fname in files:
    z = fname[fname.find(".z")+2:fname.find(".z")+5]
    if "em" in z:
        z_arr.append(10**(-1*int(z[-1])))
    else:
        z_arr.append(int(z)*1.e-3)


z_arr = np.array(z_arr)
files = np.array(files)

# sort by metallicity
inds = z_arr.argsort()
files = files[inds]
z_arr = z_arr[inds]

dat = np.genfromtxt(files[0])
age_arr = dat[:,0][:-8]

# create output array
data = np.zeros((len(age_arr),len(files)+1))
data[:,0] = age_arr

for j,fname in enumerate(files):
    dat = np.genfromtxt(fname)
    mass_frac = (dat[:,1]/dat[0,1])[:-8]
    data[:,j+1] = mass_frac


fname = out_prep+'.mass'
print("writing ",fname)
np.savetxt(fname, data,fmt='% .3f')

fname = 'zlegend.'+out_prep+'.dat'
print("writing ",fname)
np.savetxt(fname, z_arr,fmt='% .3f')

