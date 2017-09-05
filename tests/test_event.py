

def test_event():
    from zfits import FactFits

    f = FactFits('tests/resources/20160817_016.fits.fz')
    event = next(f)

    assert event is not None


def test_event_shapes():
    from zfits import FactFits

    f = FactFits('tests/resources/20160817_016.fits.fz')
    event = next(f)

    for k, v in event.items():
        print(k, v)
        if v.ndim == 1:
            assert v.shape[0] != 1
        if k in ['EventNum', 'NumBoards', 'TriggerNum', 'TriggerType']:
            assert v.ndim == 0
