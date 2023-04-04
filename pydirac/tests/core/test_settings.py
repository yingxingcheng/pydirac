from pydirac.core.settings import Settings

s = Settings()


def test_settings_add():
    s.a.a = 11
    s.a.b = 12
    s.b.a = 21
    s.b.b = 22
    print(s)
