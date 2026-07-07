'''Dummy crystal'''
import typing as t

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import pybt.tools.plotters as plotters

import xcoll as xc
from xcoll import materials

class DummyCrystal:
    '''Dummy crystal class for holding crystal parameters. Includes "dummy" tracking function'''
    def __init__(self):
        '''
        Dummy crystal constructor

        Attributes:
            dx (float): crystal width
            x0 (float): crystal center
            xp0 (float): crystal center angle
            dxp_chan (float): channeling angle
            mu_kick_chan (float): mean kick for channeling particles
            sigma_kick_chan (float): sigma kick for channeling particles
            lambda_chan (float): channeling probability
            lambda_dechan (float): dechanneling probability
            mu_kick_vr (float): mean kick for vertical radiation particles
            sigma_kick_vr (float): sigma kick for vertical radiation particles
            sigma_kick_ac (float): sigma kick for axial channeling particles
            tracker (function): function that tracks particles through the crystal

        '''
        self.dx = 1e-3
        self.dy = 0.015
        self.d0 = 0.

        self.delta0 = 0.

        self.dxp_chan = 10e-6
        self.sigma_kick_chan = 20e-6
        self.lambda_chan = 0.6
        self.lambda_dechan = 0.05  # Not implemented

        self.mu_kick_vr = -15e-6
        self.sigma_kick_vr = 20e-6

        self.sigma_kick_ac = 20e-6
        self.active_length = 0.0025
        self.bending_radius = 14.7

        self.tracker = None

    @property
    def mu_kick_chan(self):
        return self.active_length / self.bending_radius

    @mu_kick_chan.setter
    def mu_kick_chan(self, value):
        self.active_length = value * self.bending_radius

    def build_pandas_tracker(self):
        def tracker(particles):
            particles_inside_idx = np.abs(particles.x - self.d0) < self.dx
            particles_inside_size = np.sum(particles_inside_idx)
            particles_inside_px = particles.px[particles_inside_idx]
            particles_ch_idx = np.abs(particles_inside_px - self.delta0) < self.dxp_chan
            particles_ch_size = np.sum(particles_ch_idx)
            particles_vr_idx = np.abs(
                particles_inside_px - self.delta0 - self.dxp_chan - self.mu_kick_chan / 2) < self.mu_kick_chan / 2
            particles_vr_size = np.sum(particles_vr_idx)

            kicks_ch = np.random.uniform(-self.sigma_kick_chan / 2, self.sigma_kick_chan / 2, particles_ch_size)
            kicks_ch += self.mu_kick_chan * np.random.binomial(1, self.lambda_chan, particles_ch_size)
            kicks_vr = self.mu_kick_vr + np.random.uniform(-self.sigma_kick_vr / 2, self.sigma_kick_vr / 2,
                                                           particles_vr_size)
            kicks_ac = np.random.uniform(-self.sigma_kick_ac / 2, self.sigma_kick_ac / 2,
                                         particles_inside_size - particles_ch_size - particles_vr_size)

            particles.px.iloc[np.where(particles_inside_idx)[0][particles_ch_idx]] += kicks_ch
            particles.px.iloc[np.where(particles_inside_idx)[0][particles_vr_idx]] += kicks_vr
            particles.px.iloc[np.where(particles_inside_idx)[0][(~particles_ch_idx) & (~particles_vr_idx)]] += kicks_ac

        self.tracker = tracker
        return tracker

    def build_xsuite_tracker(self):
        "TODO: Check if this works on xsuite particles objects"
        def tracker(particles):
            particles_inside_idx = np.abs(particles.x - self.d0) < self.dx
            particles_inside_size = np.sum(particles_inside_idx)
            particles_inside_px = particles.px[particles_inside_idx]
            particles_ch_idx = np.abs(particles_inside_px - self.delta0) < self.dxp_chan
            particles_ch_size = np.sum(particles_ch_idx)
            particles_vr_idx = np.abs(
                particles_inside_px - self.delta0 - self.dxp_chan - self.mu_kick_chan / 2) < self.mu_kick_chan / 2
            particles_vr_size = np.sum(particles_vr_idx)

            kicks_ch = np.random.uniform(-self.sigma_kick_chan / 2, self.sigma_kick_chan / 2, particles_ch_size)
            kicks_ch += self.mu_kick_chan * np.random.binomial(1, self.lambda_chan, particles_ch_size)
            kicks_vr = self.mu_kick_vr + np.random.uniform(-self.sigma_kick_vr / 2, self.sigma_kick_vr / 2,
                                                           particles_vr_size)
            kicks_ac = np.random.uniform(-self.sigma_kick_ac / 2, self.sigma_kick_ac / 2,
                                         particles_inside_size - particles_ch_size - particles_vr_size)

            particles.px[np.where(particles_inside_idx)[0][particles_ch_idx]] += kicks_ch
            particles.px[np.where(particles_inside_idx)[0][particles_vr_idx]] += kicks_vr
            particles.px[np.where(particles_inside_idx)[0][(~particles_ch_idx) & (~particles_vr_idx)]] += kicks_ac

        self.tracker = tracker
        return tracker

    def track(self, particles, tracker: t.Literal["pandas", "xsuite"] = "pandas"):
        if tracker == "pandas":
            self.tracker = self.build_pandas_tracker()
        elif tracker == "xsuite":
            self.tracker = self.build_xsuite_tracker()
        assert self.tracker is not None, "Tracker not built"
        self.tracker(particles)

    def build_xsuite_crystal(self):
        """Builds an EverestCrystal object from dummy crystal"""
        return xc.EverestCrystal(material=materials.SiliconCrystal,
                                 #align_angle=self.xp0,  # crystal angular alignment, 0 is horizontal
                                 tilt = self.delta0,
                                 bending_angle=self.mu_kick_chan,  # bending radius [m]
                                 jaw = self.d0,  # jaw position [m] 
                                 width=self.dx*2,  # horizontal dimension of crystal [m]
                                 height=self.dy*2,  # vertical dimension of crystal [m]
                                 length=self.active_length,  # the crystal length [m]
                                 miscut=0.0,  # just set to 0 if crystal is assumed to be manufactured well
                                 lattice='strip',  # means 110
                                 #reference_center=self.x0 - self.dx,
                                 side='+',  # will only keep positive x side jaw, and will ignore negative x side jaw
                                 )