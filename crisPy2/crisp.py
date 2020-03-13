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
        plot_axis_indices : int or list, optional
            The dimensions of the cube to plot. If a list of two dimensions are given there will be an image plot with any other dimensions being used as sliders. If an int is given then that dimension will be plotted on the x-axis with the other dimensions being sliders. Default is None which implicitly passed the list [-2,-1] which are typically the x- and y-axes respectively. In the 2D case the list is parsed such that the first entry is plotted on the x-axis and the second entry is plotted on the y-axis.
        axes_coordinates : astropy.units.Quantity or list, optional
            Denotes the physical coordinates for the plots and slider axes. If None then the physical coordinates are derived from the WCS with appropriate conversions calculated. If the length equals the number of sequence dimensions i.e. the number of points, then each element describes the coordinates of the corresponding sequence dimension e.g. for n pixels in 3D space the axes_coordinates would be [(x_1,y_1,z_1), ..., (x_n,y_n,z_n)]. If the length equals the length of `plot_axis_indices`, the 0th entry describes the coordinates of the x-axis etc.; the value of each entry should be either `astropy.units.Quantity` or `numpy.ndarray` of coordinates for each pixel.
        axes_units : astropy.units.Unit or list, optional
            The units to plot on the axes. This will implicitly do the conversions to any other unit from the ones derived from the WCS. If the data is 1D then this corresponds to the units of the x-axis. If the data is 2D then this corresponds to the units of the x- and y-axis. If None the units are derived from the WCS.
        data_unit : astropy.unit.Unit, optional
            The units to plot the data in. This will implicitly do the conversions to any other unit from the ones derived from the WCS. If the data is 1D then this corresponds to the units of the y-axis. If the data is 2D then this corresponds to intensity values. If None the units are derived from the WCS.
        air : bool, optional
            This is an option to plot the spectra using the air wavelengths. This is calculated using `specutils.utils.wcs_utils.vac_to_air` with the default (Griesen 2006) method. This assumes the WCS is in terms of vacuum wavelengths. Default is False.
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