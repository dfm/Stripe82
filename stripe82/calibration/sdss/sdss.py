#!/usr/bin/env python
# encoding: utf-8
"""
Main interface to SDSS observations

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

__all__ = ['Observation','ObservationAccessError','PhotometryError']

from sdssfield import SDSSObservation
import cas

class ObservationAccessError(Exception):
    pass
class PhotometryError(Exception):
    pass

class Observation:
    """
    Wrapper class around SDSSField that makes easier interface to calibration model
    
    Parameters
    ----------
    objid : bson.ObjectID
        The field ID from the CAS MongoDB table
    
    Raises
    ------
    ObservationAccessError :
        If field can't be loaded

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey
    
    """
    def __init__(self,objid):
        self.objid = objid
        self.info  = cas.get_observation(objid)
        try:
            self.field = Observation(self.info['run'],self.info['camcol'])
            return
        except:
            raise ObservationAccessError()

    def photometry(self,ra,dec):
        """
        Return the forced photometry at RA/Dec
        
        Parameters
        ----------
        ra : float
            In degrees

        dec : float
            In degrees
        
        Returns
        -------
        counts : float
            The measured photometry in counts

        error : float
            The statistical uncertainty on counts

        Raises
        ------
        PhotometryError :
            If SDSS photometry fails

        Notes
        -----
        This is hardcoded to only use g-band for historical reasons
        
        History
        -------
        2011-06-13 - Created by Dan Foreman-Mackey
        
        """
        try:
            res,img,inv,psf = self.field.photo_at_radec(ra,dec)
        except:
            raise PhotometryError()
        return res

