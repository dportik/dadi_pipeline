import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing two population scenarios.
'''

def no_divergence(notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def anc_sym_mig(params, ns, pts):
    """
    Split with symmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
 
def anc_asym_mig(params, ns, pts):
    """
    Split with asymmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_sym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_asym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs



#######################################################################################################
#Models involving size changes

def no_mig_size(params, ns, pts):
    """
    Split with no migration, then size change with no migration.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    """
    nu1a, nu2a, nu1b, nu2b, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sym_mig_size(params, ns, pts):
    """
    Split with symmetric migration, then size change with symmetric migration.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    m: Migration rate between populations (2*Na*m)
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def asym_mig_size(params, ns, pts):
    """
    Split with different migration rates, then size change with different migration rates.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def anc_sym_mig_size(params, ns, pts):
    """
    Split with symmetrical gene flow, followed by size change with no gene flow.  

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def anc_asym_mig_size(params, ns, pts):
    """
    Split with asymmetrical gene flow, followed by size change with no gene flow.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_sym_mig_size(params, ns, pts):
    """
    Split with no gene flow, followed by size change with symmetrical gene flow.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
def sec_contact_asym_mig_size(params, ns, pts):
    """
    Split with no gene flow, followed by size change with asymmetrical gene flow.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


#######################################################################################################
#Two Epoch split with changing migration rates

def sym_mig_twoepoch(params, ns, pts):
    """
    Split into two populations, with symmetric migration. A second period of symmetric
    migration occurs, but can be a different rate. Pop size is same.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m1, m2, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m1, m21=m1)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m2, m21=m2)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def asym_mig_twoepoch(params, ns, pts):
    """
    Split into two populations, with different migration rates. A second period of asymmetric
    migration occurs, but can be at different rates. Pop size is same.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12a, m21a, m12b, m21b, T1, T2 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m12a, m21=m21a)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12b, m21=m21b)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    
    
#######################################################################################################
#Three Epoch: Divergence and Isolation, Secondary Contact, Isolation

def sec_contact_sym_mig_three_epoch(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow, then isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and third epoch.
    T3: The scaled time between the isolation and present.
    """
    nu1, nu2, m, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T3, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_asym_mig_three_epoch(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow, then isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and third epoch.
    T3: The scaled time between the isolation and present.
    """
    nu1, nu2, m12, m21, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def sec_contact_sym_mig_size_three_epoch(params, ns, pts):
    """
    Split with no gene flow, followed by size change with symmetrical gene flow, then isolation.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and isolation.
    T3: The scaled time between the isolation and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T3, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
def sec_contact_asym_mig_size_three_epoch(params, ns, pts):
    """
    Split with no gene flow, followed by size change with asymmetrical gene flow, then isolation.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and isolation.
    T3: The scaled time between the isolation and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)

    phi = Integration.two_pops(phi, xx, T3, nu1b, nu2b, m12=0, m21=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


#######################################################################################################
#Island Models - simple set
# Here populations are fractions derived from the ancestral population. The assumption is that 
# pop 2 is the 'island' population. This pop is set to s*Nref, and when these models are called 
# the upper bound on s is set to 0.5, such that the island pop can never have >50% of the 
# ancestral population (Nref = Na = 1). The values of nu1 and nu2 are therefore 1-s and s,
# (unless there is a size change event, and then it is explicit). In the absence of size change, 
# you will need to calculate these values manually, as only the s value is returned in the param list.

def vic_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. Populations are fractions of the reference 
    population, where population 2 is represented by Nref*(s), and population 1 is represented 
    by Nref*(1-s). Nref implicitly has size = 1. 
    
    params = 2
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    T, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=s, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def vic_anc_sym_mig(params, ns, pts):
    """
    Split with symmetric migration followed by isolation. Populations are fractions of the reference 
    population, where population 2 is represented by Nref*(s), and population 1 is represented 
    by Nref*(1-s). Nref implicitly has size = 1. 
    
    params = 4
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    m: Migration between pop 2 and pop 1 (2*Na*m).
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    m, T1, T2, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
   
    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=s, m12=m, m21=m)
    
    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=s, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def vic_anc_asym_mig(params, ns, pts):
    """
    Split with asymmetric migration followed by isolation. Populations are fractions of the reference 
    population, where population 2 is represented by Nref*(s), and population 1 is represented 
    by Nref*(1-s). Nref implicitly has size = 1. 
    
    params = 5
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    m12, m21, T1, T2, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
   
    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=s, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=s, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
    
def vic_sec_contact_sym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow. Populations are 
    fractions of the reference population, where population 2 is represented by Nref*(s), and 
    population 1 is represented by Nref*(1-s). Nref implicitly has size = 1. 
    
    params = 4
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    m: Migration between pop 2 and pop 1 (2*Na*m).
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    m, T1, T2, s = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=s, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=s, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def vic_sec_contact_asym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow. Populations are 
    fractions of the reference population, where population 2 is represented by Nref*(s), and 
    population 1 is represented by Nref*(1-s). Nref implicitly has size = 1. 
    
    params = 5
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    m12, m21, T1, T2, s = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=s, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=s, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_nomig(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of the reference 
    population, where population 2 is represented by Nref*(s), and population 1 is 
    represented by Nref*(1-s). Population two undergoes an exponential growth event to 
    obtain size nu2, while population one is constant. Nref implicitly has size = 1. 
	
	params = 3
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu2, T, s = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu2_func = lambda t: s * (nu2/s)**(t/T)
    
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_sym(params, ns, pts):
    """
    Split into two populations, with one migration rate. Populations are fractions of the 
    reference population, where population 2 is represented by Nref*(s), and population 1 
    is represented by Nref*(1-s). Population two undergoes an exponential growth event to 
    obtain size nu2, while population one is constant. Nref implicitly has size = 1. 
	
	params = 4
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration (2*Na*m12)
    """
    nu2, m, T, s = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu2_func = lambda t: s * (nu2/s)**(t/T)
    
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_asym(params, ns, pts):
    """
    Split into two populations, with two migration rates. Populations are fractions of the 
    reference population, where population 2 is represented by Nref*(s), and population 1 
    is represented by Nref*(1-s). Population two undergoes an exponential growth event to 
    obtain size nu2, while population one is constant. Nref implicitly has size = 1. 
	
	params = 5

    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nu2, m12, m21, T, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu2_func = lambda t: s * (nu2/s)**(t/T)

    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs




#######################################################################################################
#Island Models - discrete admixture events

# Here populations are fractions derived from the ancestral population. The assumption is that 
# pop 2 is the 'island' population. This pop is set to s*Nref, and when these models are called 
# the upper bound on s is set to 0.5, such that the island pop can never have >50% of the 
# ancestral population (Nref or Na, which = 1). The values of nu1 and nu2 are therefore 1-s and s,
# (unless there is a size change event). You will need to calculate these values manually, as only
# the s value is returned in the param list.


def vic_no_mig_admix_early(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (before drift). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
    Nref implicitly has size = 1. 
     
	params = 3
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    T: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    T, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=s, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def vic_no_mig_admix_late(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (after drift). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
    Nref implicitly has size = 1. 
    
	params = 3
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    T: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    T, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=s, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
    
def vic_two_epoch_admix(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (between two drift events). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
    Nref implicitly has size = 1. 
    
	params = 4
    
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    T1: The scaled time between the split and admixture event (in units of 2*Na generations).
    T2: The scaled time between the admixture event and present.
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    T1, T2, s, f = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=s, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=s, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_nomig_admix_early(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (before drift). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
	Population two undergoes an exponential growth event to obtain size nu2, while 
	population one is constant. Nref implicitly has size = 1. 
	
	params = 4
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    nu2, T, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    nu2_func = lambda t: s * (nu2/s)**(t/T)
    
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_nomig_admix_late(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (after drift). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
	Population two undergoes an exponential growth event to obtain size nu2, while 
	population one is constant. Admixture occurs when size nu2 is achieved. Nref implicitly has size = 1. 
	
	params = 4
	
    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    nu2, T, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu2_func = lambda t: s * (nu2/s)**(t/T)

    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def founder_nomig_admix_two_epoch(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs (between two drift events). Populations are fractions of the reference population, 
    where population 2 is represented by Nref*(s), and population 1 is represented by Nref*(1-s). 
	Population two undergoes an exponential growth event to obtain size nu2, while 
	population one is constant. Admixture occurs after size nu2 is achieved. Nref implicitly has size = 1. 
	
	params = 5

    s: Fraction of Nref that goes to pop2. Pop2 size = s. Pop 1 size = 1-s.
    nu2: Final size of pop 2, after exponential growth.
    T1: Time in the past of split (in units of 2*Na generations)
    T2: The scaled time between the admixture event and present.
    f: Fraction of updated population 2 to be derived from population 1.
    """
    nu2, T1, T2, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu2_func = lambda t: s * (nu2/s)**(t/T1)
    
    phi = Integration.two_pops(phi, xx, T1, nu1=1-s, nu2=nu2_func, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    phi = Integration.two_pops(phi, xx, T2, nu1=1-s, nu2=nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
