# beam samples

import numpy as np
import math
import scipy.optimize as opt
import scipy.special as spec
import scipy.stats as stats

from pyre.distributions import *

class RCBeam(object):
    """Some comments"""
    def __init__(self, name, ME_dict, material_dict, geo_dict):
    #def __init__(self, name, ME_flex, fc, fy, Ast, b, d, bf, hf):
        """ beam model
            :Attributs:
                ME_flex  --- model error for rc flexural strength (aci model)
                Ast --- steel area
                fy  --- yield strength
                d   --- beam depth
                bf  --- flange width
                hf  --- flange thickness
                b   --- beam width
                fc  --- compressive strength of concrete
        """
        self.name = name
        # flexural properties
        self.ME_flex = np.array(ME_dict['ME_flex'])
        self.fc = np.array(material_dict['fc'])
        self.lbda_fc = np.array(material_dict['lbda_fc'])
        self.fy = np.array(material_dict['fy'])
        self.Ast = np.array(geo_dict['Ast'])
        self.b = np.array(geo_dict['b'])
        self.d = np.array(geo_dict['d'])
        self.bf = np.array(geo_dict['bf'])
        self.hf = np.array(geo_dict['hf'])
        # shear properties
        self.ME_shear = np.array(ME_dict['ME_shear'])
        self.fyv = np.array(material_dict['fyv'])
        self.Asvt = np.array(geo_dict['Asvt'])
        self.dv = np.array(geo_dict['dv'])
        self.sv = np.array(geo_dict['sv'])

    def flexCapacity(self):
        """Flexural capacity of RC beams (aci model)
        :Returns:
            Mu --- flexural capacity (kNm)
        """
        ME = self.ME_flex
        fc = self.fc
        fy = self.fy
        Ast = self.Ast
        b = self.b
        d = self.d
        bf = self.bf
        hf = self.hf

        a = Ast * fy / (0.85* fc) / bf
        Mu = ME * Ast * fy * ( d - a/2 ) *1e-6

        beta = 0.85 - 0.05*(fc-8)/7
        beta[beta<0.65] = 0.65
        c = a/beta
        type2T = c>hf
        c[type2T] = (Ast[type2T]*fy[type2T] - 0.85*fc[type2T]*(bf[type2T]-b[type2T])*hf[type2T]) / (0.85*fc[type2T]*b[type2T]*beta[type2T])
        a[type2T] = c[type2T] * beta[type2T]
        Mu[type2T] = ME[type2T] * ( Ast[type2T] * fy[type2T] * ( d[type2T] - a[type2T]/2 ) + \
                    (0.85*fc[type2T])* (bf[type2T]-b[type2T]) * hf[type2T] * \
                    (a[type2T]/2 - hf[type2T]/2) ) *1E-6

        return Mu


    def shearCapacity(self):

        fc = self.fc
        b = self.b
        dv = self.dv
        Asvt = self.Asvt
        fyv = self.fyv
        sv = self.sv
        ME_shear = self.ME_shear

        Vc = 0.17 * self.lbda_fc * np.sqrt(fc) * b * dv
        Vs = Asvt * fyv * dv / sv
        Vs_max = 0.66 * np.sqrt(fc) * b * dv
        over_reinf = Vs > Vs_max
        Vs[over_reinf] = Vs_max[over_reinf]

        Vu = ME_shear*(Vc + Vs)*1e-3    # [kN]

        return Vu
