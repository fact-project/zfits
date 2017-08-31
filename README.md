# zfits


from zfits import FactFits
with FactFits('20160401_153.fits.fz') as f:
    for event in f:
        print(event)

# event is a dict{string: np.1d-array}


# Build with:

    export CFLAGS='-std=c++11' && pip install -e .
