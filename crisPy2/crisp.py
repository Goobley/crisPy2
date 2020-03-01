import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from specutils.utils.wcs_utils import vac_to_air
from ndcube import NDCube, NDCubeSequence
from astropy.wcs import WCS
import astropy.units as u
import html

class CRISPSpectro(NDCube):
    '''
    This is the class for a singular CRISP cube described by one WCS.
    '''

    def __init__(self, data, wcs=None, meta=None, uncertainty=None, mask=None, unit=None, extra_coords=None, copy=False, missing_axes=None):
        super().__init__(data, wcs, meta=meta, uncertainty=uncertainty, mask=mask, unit=unit, extra_coords=extra_coords, copy=copy, missing_axes=missing_axes)
        self.aa = html.unescape("&#8491;")
        self.l = html.unescape("&lambda;")

    def __str__(self):
        if self.meta is None:
            return "You know more about this data than I do."
        else:
            time = self.meta.get("DATE-AVG")[-12:]
            date = self.meta.get("DATE-AVG")[:-13]
            cl = str(np.round(self.meta.get("TWAVE1")))
            wwidth = self.meta.get("WWIDTH1")
            shape = str([self.meta.get("NAXIS3"), self.meta.get("NAXIS2"), self.meta.get("NAXIS1")])
            el = self.meta.get("WDESC1")
            pointing_x = str(self.meta.get("CRVAL1"))
            pointing_y = str(self.meta.get("CRVAL2"))

            return f"""CRISP Imaging Spectroscopy
            -------------------------
            {date} {time}
            
            Observed: {el}
            Centre wavelength: {cl}
            Wavelengths sampled: {wwidth}
            Pointing: ({pointing_x}, {pointing_y})
            Shape: {shape}"""

    def __repr__(self):
        return self.__str__()

    def plot(self, axes=None, plot_axis_indices=None, axes_coordinates=None, axes_units=None, data_unit=None, **kwargs):
        if len(self.dimensions) == 1:
            super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=axes_coordinates, data_unit=data_unit, ylabel="Intensity [DNs]", xlabel=f"{self.l} [{self.aa}]", **kwargs)
        else:
            ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=axes_coordinates, data_unit=data_unit, **kwargs)
            ax.set_ylabel("Helioprojective Latitude [arcsec]")
            ax.set_xlabel("Helioprojective Longitude [arcsec]")
            ax.grid()
            ax.colorbar()