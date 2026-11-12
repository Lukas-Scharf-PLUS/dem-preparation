"""
Write the final DEM in ISCE2 format: ``<name>.dem.wgs84`` (raw float32),
``<name>.dem.wgs84.xml`` (ISCE header) and ``<name>.dem.wgs84.vrt`` (GDAL).

This is the format ``stackSentinel.py -d`` expects. ISCE2 is not needed to
write it; the XML has the same content as the one ISCE2's ``renderHdr()``
writes for a DEM image.

Pixel convention: ISCE2's ``startingValue`` is the *center* of the first
pixel, GDAL's geotransform the *corner*. Half a pixel is added for the XML.
(``gdal2isce_xml.py`` and ISCE2's own VRT writer skip this and shift the DEM
by half a pixel; the VRT written here uses the correct corner.)

ISCE2 has no no-data value for DEMs. No-data pixels are written as 0, as in
ISCE2's own DEM stitcher (there should be none if the coverage check passed).
"""

import logging

import numpy as np
from osgeo import gdal

logger = logging.getLogger(__name__)
gdal.UseExceptions()

ISCE_SUFFIX = ".dem.wgs84"
ISCE_NODATA_HEIGHT = 0.0
ROWS_PER_BLOCK = 2048

XML_TEMPLATE = """<imageFile>
    <property name="access_mode"><value>read</value></property>
    <property name="byte_order"><value>l</value></property>
    <component name="coordinate1">
        <factorymodule>isceobj.Image</factorymodule>
        <factoryname>createCoordinate</factoryname>
        <doc>First coordinate of a 2D image (width).</doc>
        <property name="delta"><value>{delta_lon!r}</value></property>
        <property name="endingvalue"><value>{end_lon!r}</value></property>
        <property name="family"><value>imagecoordinate</value></property>
        <property name="name"><value>imagecoordinate_name</value></property>
        <property name="size"><value>{width}</value></property>
        <property name="startingvalue"><value>{first_lon!r}</value></property>
    </component>
    <component name="coordinate2">
        <factorymodule>isceobj.Image</factorymodule>
        <factoryname>createCoordinate</factoryname>
        <doc>Second coordinate of a 2D image (length).</doc>
        <property name="delta"><value>{delta_lat!r}</value></property>
        <property name="endingvalue"><value>{end_lat!r}</value></property>
        <property name="family"><value>imagecoordinate</value></property>
        <property name="name"><value>imagecoordinate_name</value></property>
        <property name="size"><value>{length}</value></property>
        <property name="startingvalue"><value>{first_lat!r}</value></property>
    </component>
    <property name="data_type"><value>FLOAT</value></property>
    <property name="family"><value>demimage</value></property>
    <property name="file_name"><value>{file_name}</value></property>
    <property name="image_type"><value>dem</value></property>
    <property name="length"><value>{length}</value></property>
    <property name="name"><value>demimage_name</value></property>
    <property name="number_bands"><value>1</value></property>
    <property name="reference"><value>WGS84</value></property>
    <property name="scheme"><value>BIP</value></property>
    <property name="width"><value>{width}</value></property>
    <property name="xmax"><value>{end_lon!r}</value></property>
    <property name="xmin"><value>{first_lon!r}</value></property>
</imageFile>
"""

VRT_TEMPLATE = """<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <SRS>EPSG:4326</SRS>
    <GeoTransform>{west!r}, {delta_lon!r}, 0.0, {north!r}, 0.0, {delta_lat!r}</GeoTransform>
    <VRTRasterBand dataType="Float32" band="1" subClass="VRTRawRasterBand">
        <SourceFilename relativeToVRT="1">{file_name}</SourceFilename>
        <ByteOrder>LSB</ByteOrder>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>4</PixelOffset>
        <LineOffset>{line_offset}</LineOffset>
    </VRTRasterBand>
</VRTDataset>
"""


def isce_header_values(geotransform, width, length):
    """Return the ISCE header values (pixel-center convention) for a GDAL geotransform."""
    west, delta_lon, _, north, _, delta_lat = geotransform
    first_lon = west + 0.5 * delta_lon
    first_lat = north + 0.5 * delta_lat
    return {
        "first_lon": first_lon,
        "first_lat": first_lat,
        "delta_lon": delta_lon,
        "delta_lat": delta_lat,
        # Same definition as ISCE2: start + size * delta.
        "end_lon": first_lon + width * delta_lon,
        "end_lat": first_lat + length * delta_lat,
        "width": width,
        "length": length,
    }


def write_raw_heights(dataset, out_file):
    """Copy band 1 of ``dataset`` to a raw little-endian float32 file, block by block.

    Returns
    -------
    int
        Number of no-data pixels that were written as 0.
    """
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    replaced = 0
    with open(out_file, "wb") as stream:
        for first_row in range(0, dataset.RasterYSize, ROWS_PER_BLOCK):
            rows = min(ROWS_PER_BLOCK, dataset.RasterYSize - first_row)
            block = band.ReadAsArray(0, first_row, dataset.RasterXSize, rows).astype("<f4")
            invalid = ~np.isfinite(block)
            if nodata is not None:
                invalid |= block == nodata
            replaced += int(invalid.sum())
            block[invalid] = ISCE_NODATA_HEIGHT
            block.tofile(stream)
    return replaced


def export_isce_dem(source_file, out_base):
    """Write ``source_file`` (lat/lon GDAL raster) as an ISCE2 DEM.

    Parameters
    ----------
    source_file : Path
        Final DEM (GeoTIFF or VRT), EPSG:4326, ellipsoid heights.
    out_base : Path
        Output path without suffix; ``.dem.wgs84`` (+ ``.xml``, ``.vrt``) is appended.

    Returns
    -------
    Path
        The raw ``.dem.wgs84`` file.
    """
    out_file = out_base.with_name(out_base.name + ISCE_SUFFIX)
    dataset = gdal.Open(str(source_file))
    width, length = dataset.RasterXSize, dataset.RasterYSize
    geotransform = dataset.GetGeoTransform()

    logger.info("Writing ISCE2 DEM %s (%d x %d pixels, %.1f GB)", out_file, width, length,
                width * length * 4 / 1e9)
    replaced = write_raw_heights(dataset, out_file)
    if replaced:
        logger.warning("%d no-data pixels written as %g m in the ISCE DEM.", replaced, ISCE_NODATA_HEIGHT)

    values = isce_header_values(geotransform, width, length)
    # Absolute path, like ISCE2 writes it; stackSentinel.py reads the XML next to the file anyway.
    xml_text = XML_TEMPLATE.format(file_name=out_file.resolve(), **values)
    out_file.with_name(out_file.name + ".xml").write_text(xml_text, encoding="utf-8")
    vrt_text = VRT_TEMPLATE.format(
        width=width, length=length, west=geotransform[0], delta_lon=geotransform[1],
        north=geotransform[3], delta_lat=geotransform[5], file_name=out_file.name,
        line_offset=4 * width,
    )
    out_file.with_name(out_file.name + ".vrt").write_text(vrt_text, encoding="utf-8")
    return out_file
