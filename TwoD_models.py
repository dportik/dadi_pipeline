import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing two population scenarios.
First three models taken from Demographics2D.py in dadi package.
Several models follow examples from Rougemont et al 2016 (doi:10.1111/mec.13664).
Instantaneous size change is added to several of these.
'''

def snm(notused, ns, pts):
    """
    ns = (n1,n2)

    Standard neutral model, populations never diverge.
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def bottlegrowth(params, ns, pts):
    """
    params = (nuB,nuF,T)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth with no population
    split.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T = params
    return bottlegrowth_split_mig((nuB,nuF,0,T,0), ns, pts)

def bottlegrowth_split(params, ns, pts):
    """
    params = (nuB,nuF,T,Ts)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth then split.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T,Ts = params
    return bottlegrowth_split_mig((nuB,nuF,0,T,Ts), ns, pts)

def bottlegrowth_split_mig(params, ns, pts):
    """
    params = (nuB,nuF,m,T,Ts)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth then split with
    migration.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    m: Migration rate between the two populations (2*Na*m).
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,m,T,Ts = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = Integration.one_pop(phi, xx, T-Ts, nu_func)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu0 = nu_func(T-Ts)
    nu_func = lambda t: nu0*numpy.exp(numpy.log(nuF/nu0) * t/Ts)
    phi = Integration.two_pops(phi, xx, Ts, nu_func, nu_func, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_mig(params, ns, pts):
    """
    params = (nu1,nu2,T,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_asym_mig(params, ns, pts):
    """
    params = (nu1,nu2,T,m12,m21)
    ns = (n1,n2)

    Split into two populations of specifed size, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
	"""
    nu1, nu2, T, m12, m21 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def split_no_mig(params, ns, pts):
    """
    params = (nu1,nu2,T)
    ns = (n1,n2)

    Split into two populations of specifed size, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

 
def split_ancient_asymig(params, ns, pts):

    """
    Model with split and migration, period of no size change and no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1, nu2, m12, m21, T, Tam = params

    xx = Numerics.default_grid(pts)


    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)

    
    phi = Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_ancient_symig(params, ns, pts):

    """
    Model with split and migration, period of no size change and no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1, nu2, m, T, Tam = params

    xx = Numerics.default_grid(pts)


    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    
    phi = Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_secondary_contact_asym(params, ns, pts):

    """
    Model with split and no gene flow, followed by period of no size change and asymmetrical gene flow

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1, nu2, m12, m21, T, Tsc = params

    xx = Numerics.default_grid(pts)

    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)


    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_secondary_contact_sym(params, ns, pts):

    """
    Model with split and no gene flow, followed by period of no size change and symmetrical gene flow

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1, nu2, m, T, Tsc = params

    xx = Numerics.default_grid(pts)

    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m, m21=m)


    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


#######################################################################################################
#Models involving instantaneous size changes

def split_no_mig_size(params, ns, pts):
    """
    params = (nu1b,nu2b,nu1r,nu2r,T1,T2)
    ns = (n1,n2)

    Split into two populations of specifed size, with no migration, period of size change and no migration.

    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1r: Size of population 1 after time interval.
    nu2r: Size of population 2 after time interval.
    T2: Time of population size change.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1b,nu2b,nu1r,nu2r,T1,T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1b, nu2b, m12=0, m21=0)
    
    phi = Integration.two_pops(phi, xx, T2, nu1r, nu2r, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

 
def split_ancient_asymmig_size(params, ns, pts):

    """
    Model with split and no gene flow, followed by period of size change and asymmetrical gene flow.

    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1r: Size of population 1 after time interval.
    nu2r: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1b,nu2b,nu1r,nu2r,T1,T2,m12,m21 = params

    xx = Numerics.default_grid(pts)


    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1b, nu2b, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1r, nu2r, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_ancient_symmig_size(params, ns, pts):

    """
    Model with split and no gene flow, followed by period of size change and symmetrical gene flow.

    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1r: Size of population 1 after time interval.
    nu2r: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m: Migration between pop 2 and pop 1.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1b,nu2b,nu1r,nu2r,T1,T2,m = params

    xx = Numerics.default_grid(pts)


    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1b, nu2b, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1r, nu2r, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_secondary_contact_asym_size(params, ns, pts):

    """
    Model with split, complete isolation, followed by secondary contact with asymmetrical gene flow

    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1r: Size of population 1 after time interval.
    nu2r: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1b,nu2b,nu1r,nu2r,T1,T2,m12,m21 = params

    xx = Numerics.default_grid(pts)

    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1r, nu2r, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_secondary_contact_sym_size(params, ns, pts):

    """
    Model with split, complete isolation, followed by secondary contact with asymmetrical gene flow

    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1r: Size of population 1 after time interval.
    nu2r: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m: Migration between pop 2 and pop 1.
    ns: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1b,nu2b,nu1r,nu2r,T1,T2,m = params

    xx = Numerics.default_grid(pts)

    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1r, nu2r, m12=m, m21=m)


    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
