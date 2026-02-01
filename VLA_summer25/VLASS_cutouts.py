import numpy as np
import os
import sys
import argparse
import glob
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from urllib.request import urlopen
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.time import Time
import ssl
from ztfquery.utils import stamps
import requests
from PIL import Image
import io
import csv
import logging
import urllib.request


ssl._create_default_https_context = ssl._create_unverified_context

fname = "VLASS_dyn_summary.php"
url = 'https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php'
output_file = 'CSV'

urllib.request.urlretrieve(url, output_file)
print(f'File downloaded to: {output_file}')


def get_tiles():
    """ Get tiles
    I ran wget https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php
    """

    inputf = open(fname, "r")
    lines = inputf.readlines()
    inputf.close()

    header = list(filter(None, lines[0].split("  ")))

    header = np.array([val.strip() for val in header])

    names = []
    dec_min = []
    dec_max = []
    ra_min = []
    ra_max = []
    obsdate = []
    epoch = []

    for line in lines[3:]:
        dat = list(filter(None, line.split("  ")))
        dat = np.array([val.strip() for val in dat])
        names.append(dat[0])
        dec_min.append(float(dat[1]))
        dec_max.append(float(dat[2]))
        ra_min.append(float(dat[3]))
        ra_max.append(float(dat[4]))
        obsdate.append(dat[6])
        epoch.append(dat[5])

    names = np.array(names)
    dec_min = np.array(dec_min)
    dec_max = np.array(dec_max)
    ra_min = np.array(ra_min)
    ra_max = np.array(ra_max)
    obsdate = np.array(obsdate)
    epoch = np.array(epoch)

    return names, dec_min, dec_max, ra_min, ra_max, epoch, obsdate


def search_tiles(tiles, c):
    """ Now that you've processed the file, search for the given RA and Dec

    Parameters
    ----------
    c: SkyCoord object
    """
    ra_h = c.ra.hour
    dec_d = c.dec.deg
    names, dec_min, dec_max, ra_min, ra_max, epochs, obsdate = tiles
    has_dec = np.logical_and(dec_d > dec_min, dec_d < dec_max)
    has_ra = np.logical_and(ra_h > ra_min, ra_h < ra_max)
    in_tile = np.logical_and(has_ra, has_dec)
    name = names[in_tile]
    epoch = epochs[in_tile]
    date = obsdate[in_tile]
    if len(name) == 0:
        print("Sorry, no tile found.")
        return None, None, None
    else:
        return name, epoch, date


def get_subtiles(tilename, epoch):
    """ For a given tile name, get the filenames in the VLASS directory.
    Parse those filenames and return a list of subtile RA and Dec.
    RA and Dec returned as a SkyCoord object
    """
    if epoch == 'VLASS1.2':
        epoch = 'VLASS1.2v2'
    elif epoch == 'VLASS1.1':
        epoch = 'VLASS1.1v2'
    url_full = 'https://archive-new.nrao.edu/vlass/quicklook/%s/%s/' % (epoch, tilename)
    print(url_full)
    urlpath = urlopen(url_full)

    string = (urlpath.read().decode('utf-8')).split("\n")

    vals = np.array([val.strip() for val in string])

    keep_link = np.array(["href" in val.strip() for val in string])

    keep_name = np.array([tilename in val.strip() for val in string])

    string_keep = vals[np.logical_and(keep_link, keep_name)]

    fname = np.array([val.split("\"")[7] for val in string_keep])

    pos_raw = np.array([val.split(".")[4] for val in fname])

    ra = []
    dec = []

    for ii in range(len(pos_raw)):
        rah = pos_raw[ii][1:3]
        if rah == "24":
            rah = "00"
        ram = pos_raw[ii][3:5]
        ras = pos_raw[ii][5:7]
        decd = pos_raw[ii][7:10]
        decm = pos_raw[ii][10:12]
        decs = pos_raw[ii][12:]
        hms = "%sh%sm%ss" % (rah, ram, ras)
        ra.append(hms)
        dms = "%sd%sm%ss" % (decd, decm, decs)
        dec.append(dms)
    ra = np.array(ra)
    dec = np.array(dec)
    c_tiles = SkyCoord(ra, dec, frame='icrs')
    return fname, c_tiles

save_directory = "/Users/Djslime07/VLA_Summer25/newer_images"
def get_cutout(imname, name, c, epoch):
    print("Generating cutout")
    # Position of source
    ra_deg = c.ra.deg
    dec_deg = c.dec.deg

    print("Cutout centered at position %s, %s" % (ra_deg, dec_deg))

    # Open image and establish coordinate system
    try:
        with pyfits.open(imname, ignore_missing_simple=True) as hdulist:
            im = hdulist[0].data[0, 0]
            if im.size == 0:
                print("Error: Image data is empty.")
                return None

            w = WCS(hdulist[0].header)
    except Exception as e:
        print(f"Error reading image data: {e}")
        return None

    # Find the source position in pixels
    src_pix = w.wcs_world2pix([[ra_deg, dec_deg, 0, 0]], 0)
    x = src_pix[0, 0]
    y = src_pix[0, 1]

    # Check if the source is actually in the image
    pix1 = hdulist[0].header['CRPIX1']
    pix2 = hdulist[0].header['CRPIX2']
    badx = np.logical_or(x < 0, x > 2 * pix1)
    bady = np.logical_or(y < 0, y > 2 * pix2)
    if np.logical_and(badx, bady):
        print("Tile has not been imaged at the position of the source")
        return None
    else:
        # Set the dimensions of the image
        image_dim_arcsec = 12
        delt1 = hdulist[0].header['CDELT1']
        delt2 = hdulist[0].header['CDELT2']
        cutout_size = image_dim_arcsec / 3600  # in degrees
        dside1 = -cutout_size / 2. / delt1
        dside2 = cutout_size / 2. / delt2

        vmin = -1e-4
        vmax = 1e-3

        im_plot_raw = im[int(y - dside1):int(y + dside1), int(x - dside2):int(x + dside2)]
        im_plot = np.ma.masked_invalid(im_plot_raw)

        # 3-sigma clipping
        rms_temp = np.ma.std(im_plot)
        keep = np.ma.abs(im_plot) <= 3 * rms_temp
        rms = np.ma.std(im_plot[keep])

        if im_plot.flatten().size == 0:
            print("Tile has not been imaged at the position of the source")
            return None
        else:
            peak_flux = np.ma.max(im_plot.flatten())


        save_path = os.path.join(save_directory, f"{name}_{epoch}.png")

        fig, ax = plt.subplots(figsize=(6,6))
        ax.imshow(np.flipud(im_plot),
            extent=[-0.5 * cutout_size * 3600., 0.5 * cutout_size * 3600.,
                    -0.5 * cutout_size * 3600., 0.5 * cutout_size * 3600],
            vmin=vmin, vmax=vmax, cmap='YlOrRd')

        peakstr = "Peak Flux %s mJy" % (np.round(peak_flux * 1e3, 3))
        rmsstr = "RMS Flux %s mJy" % (np.round(rms * 1e3, 3))
        obsdate = hdulist[0].header['DATE-OBS']

        plt.title(f"{name} {epoch}\n{peakstr}, {rmsstr} \n {obsdate}",)
        plt.xlabel("Offset in RA (arcsec)")
        plt.ylabel("Offset in Dec (arcsec)")

        plt.savefig(save_path)
        plt.close()
        print("PNG Downloaded Successfully")

        return peak_flux, rms


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
csv_file_path = "/Users/Djslime07/VLA_Summer25/temp2.csv"
ddir = "/Users/Djslime07/VLA_Summer25/z"


def run_search(name, c, date=None):
    """
    Searches the VLASS catalog for a source

    Parameters
    ----------
    name: name of the sources
    c: coordinates as SkyCoord object
    date: date in astropy Time format
    """
    print("Running for %s" % name)
    print("Coordinates %s" % c)
    print("Date: %s" % date)

    # Find the VLASS tile(s)
    tiles = get_tiles()
    tilenames, epochs, obsdates = search_tiles(tiles, c)

    past_epochs = ["VLASS1.1v2", "VLASS1.2v2", "VLASS2.1", "VLASS2.2", "VLASS3.1", "VLASS3.2"]

    if tilenames[0] is None:
        print("There is no VLASS tile at this location")

    else:
        for ii, tilename in enumerate(tilenames):
            print()
            print("Looking for tile observation for %s" % tilename)
            epoch = epochs[ii]
            obsdate = obsdates[ii]

            valid_obsdate = True
            try:
                obsdate_mjd = Time(obsdate, format='iso').mjd
            except ValueError:
                print(f"Could not parse observation date '{obsdate}' for tile {tilename} in epoch {epoch}")
            valid_obsdate = False

            # Adjust name so it works with the version 2 ones for 1.1 and 1.2
            if epoch == 'VLASS1.2':
                epoch = 'VLASS1.2v2'
            elif epoch == 'VLASS1.1':
                epoch = 'VLASS1.1v2'

            if epoch not in past_epochs:
                    # Make list of observed tiles
                    url_full = 'https://archive-new.nrao.edu/vlass/quicklook/%s/' % (epoch)
                    urlpath = urlopen(url_full)
                    # Get site HTML coding
                    string = (urlpath.read().decode('utf-8')).split("\n")
                    # clean the HTML elements of trailing and leading whitespace
                    vals = np.array([val.strip() for val in string])
                    # Make list of useful html elements
                    files = np.array(['alt="[DIR]"' in val.strip() for val in string])
                    useful = vals[files]
                    # Splice out the name from the link
                    obsname = np.array([val.split("\"")[7] for val in useful])
                    observed_current_epoch = np.char.replace(obsname, '/', '')

                    # Check if tile has been observed yet for the current epoch
                    if epoch not in observed_current_epoch:
                        print("Sorry, tile will be observed later in this epoch")
            else:
                print("Tile Found:")
                print(tilename, epoch)
                subtiles, c_tiles = get_subtiles(tilename, epoch)
                # Find angular separation from the tiles to the location
                dist = c.separation(c_tiles)
                # Find tile with the smallest separation
                subtile = subtiles[np.argmin(dist)]
                url_get = "https://archive-new.nrao.edu/vlass/quicklook/%s/%s/%s" % (
                    epoch, tilename, subtile)
                imname = "%s.I.iter1.image.pbcor.tt0.subim.fits" % subtile[0:-1]
                fname = url_get + imname
                print(fname)
                if len(glob.glob(imname)) == 0:
                    cmd = "curl -O %s" % fname
                    print(cmd)
                    os.system(cmd)
                # Get image cutout and save FITS data as png
                out = get_cutout(imname, name, c, epoch)
                if out is not None:
                    peak, rms = out
                    output_file = "/Users/Djslime07/VLA_Summer25/fluxes_and_rms2.csv"
                    with open(output_file, 'a') as f:
                        print(f"{name}_{epoch}.png has a peak flux of {np.round(peak * 1e3, 3)} mJy and an RMS of {np.round(rms * 1e3, 3)} mJy", file=f)

                if valid_obsdate:
                    with open(output_file, 'a') as f:
                        print(f"{name}_{epoch}.png observed on {obsdate_mjd}", file=f)

                # Remove FITS file
                #os.remove(imname)

def plot_ls_cutout(ddir, name, ra_str, dec_str, outputf):
    """ Plot cutout from Legacy Survey """
    # Create a SkyCoord object from ra_str and dec_str
    coord = SkyCoord(ra_str, dec_str, unit=("hour", "degree"))

    fname = os.path.join(ddir, f"{name}_LegSurvey.png")

    if not os.path.isfile(fname):
        url = f"http://legacysurvey.org/viewer/cutout.jpg?ra={coord.ra.deg}&dec={coord.dec.deg}&layer=ls-dr9&pixscale=0.27&bands=grz"
        plt.figure(figsize=(2.1, 2.1), dpi=120)
        try:
            r = requests.get(url)
            img = Image.open(io.BytesIO(r.content))
            plt.imshow(img)
            plt.title("LegSurv DR9", fontsize=12)
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(fname, bbox_inches="tight")

            decsign = '+' if coord.dec.deg >= 0 else '-'
            lslinkstr = f"http://legacysurvey.org/viewer?ra={coord.ra.deg:.6f}&dec={decsign}{abs(coord.dec.deg):.6f}&zoom=16&layer=dr9"
            outputf.write(f"<a href='{lslinkstr}'>")
            outputf.write(f'<img src="{name}_LegSurvey.png" height="200">')
            outputf.write("</a>")
            outputf.write('<br>')
        except Exception as e:
            # Not in footprint or another error
            print(f"Error: {e}")
            return None
        finally:
            plt.close()

    return fname

def plot_ps1_cutout(ddir, name, ra_str, dec_str, outputf):
    """ Plot cutout from PanSTARRS """
    # Create a SkyCoord object from ra_str and dec_str
    coord = SkyCoord(ra_str, dec_str, unit=("hour", "degree"))

    fname = os.path.join(ddir, f"{name}_ps1.png")

    if not os.path.isfile(fname):
        img = stamps.get_ps_stamp(coord.ra.deg, coord.dec.deg, size=240, color=["y", "g", "i"])
        plt.figure(figsize=(2.1, 2.1), dpi=120)
        plt.imshow(np.asarray(img))
        plt.title("PS1 (y/g/i)", fontsize=12)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fname, bbox_inches="tight")

        decsign = '+' if coord.dec.deg >= 0 else '-'
        pslinkstr = f"https://ps1images.stsci.edu/cgi-bin/ps1cutouts?ra={coord.ra.deg:.6f}&dec={decsign}{abs(coord.dec.deg):.6f}&size=240&format=jpeg&filters=ygi"
        outputf.write(f"<a href='{pslinkstr}'>")
        outputf.write(f'<img src="{name}_ps1.png" height="200">')
        outputf.write("</a>")
        outputf.write('<br>')
        plt.close()

    return fname
def process_object(row):
    """Process a single object given a row from the CSV."""
    name = row['name']
    ra = row['ra']
    dec = row['dec']
    Obj = SkyCoord(ra, dec, unit=("hourangle", "deg"))

    try:
        logging.info(f"Processing object {name} at RA: {ra}, Dec: {dec}")
        run_search(name, Obj)

        with open("output.html", "a") as outputf:
            plot_ls_cutout(ddir, name, ra, dec, outputf)
            plot_ps1_cutout(ddir, name, ra, dec, outputf)

    except Exception as e:
        logging.error(f"Failed to process object {name}: {e}")

def process_csv(csv_file_path, start_line=0, end_line=None):
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for idx, row in enumerate(reader):
            if idx < start_line:
                continue
            if end_line is not None and idx > end_line:
                break
            process_object(row)

start_line = 0  # If you want line x, then do (x-2) for the actual line
end_line = 17   # The last line to process -2

process_csv(csv_file_path, start_line, end_line)

logging.info("Processing complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= \
                                         '''
                                         Searches VLASS for a source.
                                         User needs to supply name, RA (in decimal degrees),
                                         Dec (in decimal degrees), and (optionally) date (in mjd).
                                         If there is a date, then will only return VLASS images taken after that date
                                         (useful for transients with known explosion dates).
                                 
                                         Usage: vlass_search.py <Name> <RA [deg]> <Dec [deg]> <(optional) Date [mjd]>
                                         ''', formatter_class=argparse.RawTextHelpFormatter)

    if len(sys.argv) < 3:
        print("Usage: vlass_search.py <Name> <RA [deg]> <Dec [deg]> <(optional) Date [astropy Time]>")
        sys.exit()

    name = str(sys.argv[1])
    ra = float(sys.argv[2])
    dec = float(sys.argv[3])
    c = SkyCoord(ra, dec, unit='deg')

    if glob.glob("/Users/annaho/Dropbox/astro/tools/Query_VLASS/VLASS_dyn_summary.php"):
        pass
    else:
        print("Tile summary file is not here. Download it using wget:\
               wget https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php")

    if (len(sys.argv) > 4):
        date = Time(float(sys.argv[4]), format='mjd')
        print('Searching for observations after %s' % date)
        run_search(name, c, date)
    else:
        print('Searching all obs dates')
        run_search(name, c)