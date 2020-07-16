# Hello `cathy`!

This is a command line application for manipulating catheter data.

## Commands

Documentation about specific parameters and options should be available by examining each specific command's generated help documentation.

### `cathy info`:
Get information about a .projection file, or a directory containing related projections files.

### `cathy peek`:
Visualize catheter projections with optional coil location data.

### `cathy localize`:
Apply localization algorithms to process projections files into catheter coordinates files.

### `cathy apply-ukf`:
Smooth and fit catheter coordinate data using an unscented Kalman filter developed for this purpose.

### `cathy scatter`
Visualize catheter coordinates on a dicom image.

## Resources

### Click
Click is the library used to generate the command line interface for cathy.
It generates a CLI from a tree of functions with appropriate decorators. Documentation and extensive examples can be found here:
https://click.palletsprojects.com/en/7.x/.