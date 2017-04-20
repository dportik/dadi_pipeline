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
    Model with split and symmetric migration followed by isolation.

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
    Model with split and asymmetric migration followed by isolation.

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
    Model with split and no gene flow, followed by period of symmetrical gene flow.

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
    Model with split and no gene flow, followed by period of asymmetrical gene flow.

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
    Split into two populations with no migration, then size change.

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
    Split into two populations with symmetric migration, size change.

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
    Split into two populations with different migration rates, size change.

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
    Model with split and no gene flow, followed by symmetrical gene flow, size change.

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
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

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
    Model with split in isolation, followed by secondary contact with asymmetrical gene flow, size change.

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
    Model with split in isolation, followed by secondary contact with asymmetrical gene flow, size change.

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
