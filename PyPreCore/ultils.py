from numpy import sqrt, sum

def magnitude(r, keepdims = False):
    """ 
    scalar magnitude of a vector
    from num astro
    """
    return sqrt(sum(r**2,1,keepdims=keepdims))

    