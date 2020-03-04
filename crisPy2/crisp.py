import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from specutils.utils.wcs_utils import vac_to_air, air_to_vac
from ndcube import NDCube, NDCubeSequence
from astropy.wcs import WCS
import astropy.units as u
import html

class CRISPSpectro(NDCube):
    '''
    This is the class for a singular CRISP cube described by one WCS.

    Parameters
    ----------
    data : numpy.ndarray
        The data to be stored in the CRISPSpectro object.
    wcs : astropy.wcs.WCS, optional
        The world coordinate system to describe the observation by which can be custom or derived from the observation's metadata.
    meta : dict, optional
        The metadata for the observation.
    uncertainty : numpy.ndarray, optional
        Uncertainty array with same dimensions as data array containing the uncertainty in the observations.
    mask : numpy.ndarray, optional
        Mask for the observation following the `numpy convention <https://docs.scipy.org/doc/numpy/reference/maskedarray.html>`_.
    unit : Unit-like, optional
        Unit for the data in the observation.
    extra_coords : list, optional
        One can add an extra axis/extra axes if required by the data but not specified by the world coordinate system. `À la <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_.
    copy : bool, optional
        Whether or not to create copy of everything e.g. load it into memory with its own reference. Default is False.
    missing_axes : list, optional
        This tells us about missing axes in the data if there is any. True corresponds to the axis being missing, False corresponds to the axis being a data axis. Default is None.
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

    def plot(self, axes=None, plot_axis_indices=None, axes_coordinates=None, axes_units=None, data_unit=None, air=False, **kwargs):
        '''
        This function plots 1D or 2D slices of the data with correct axes etc.

        Parameters
        ----------
        axes : astropy.visualization.wcsaxes.core.WCSAxes, optional
            The axes to plot onto. Default is None i.e. use axes defined by this object. This kwarg is typically useful for plotting multiple plots on one axes.
        plot_axis_indices : 
        '''

        #TODO: Add set_ylabel and set_xlabel after plot for 1D to see if that fixes things.
        if (len(self.dimensions) == 1):
            if not air:
                if axes_units is None:
                    ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=axes_coordinates, data_unit=data_unit, axes_units=u.Angstrom, ylabel="Intensity [DNs]", xlabel=f"Vacuum {self.l} [{self.aa}]", fmt="o", **kwargs)
                else:
                    ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=axes_coordinates, data_unit=data_unit, axes_units=axes_units, ylabel="Intensity [DNs]", xlabel=f"Vacuum {self.l} [{axes_units}]", fmt="o", **kwargs)
            else:
                air_wvls = vac_to_air(self.axis_world_coords(0))
                if axes_units is None:
                    ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=air_wvls, data_unit=data_unit, axes_units=u.Angstrom, ylabel="Intensity [DNs]", xlabel=f"Air {self.l} [{self.aa}]", fmt="o", **kwargs)
                else:
                    ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=air_wvls, data_unit=data_unit, axes_units=axes_units, ylabel="Intensity [DNs]", xlabel=f"Air {self.l} [{self.aa}]", fmt="o", **kwargs)
            ax.set_title(self.meta.get("WDESC1")+self.aa)
        else:
            ax = super().plot(axes=axes, plot_axis_indices=plot_axis_indices, axes_coordinates=axes_coordinates, data_unit=data_unit, axes_units=axes_units, **kwargs)
            ax.set_ylabel("Helioprojective Latitude [arcsec]")
            ax.set_xlabel("Helioprojective Longitude [arcsec]")
            ax.grid()
            plt.colorbar(ax.get_images()[0], ax=ax, label="Intensity [DNs]")