version = "25.2.1"

def get_version():
    return version

def safe_int(s):
    # Remove trailing non-digits
    while s and not s.isdigit():
        s = s[:-1]
    try:
        return int(s)
    except:
        return 9999

def parse_version(v):
    parts = v.split('.')
    if len(parts) == 3:
        a,b,c = parts
        return (safe_int(a), safe_int(b), safe_int(c))
    elif len(parts) == 2:
        a, b = parts
        return (safe_int(a), safe_int(b), 0)
    else:
        return (9999,9999,9999)

GIT_VERSION = 'GIT_DESCRIBE'

if __name__ == '__main__':
    print(get_version(), end='')
