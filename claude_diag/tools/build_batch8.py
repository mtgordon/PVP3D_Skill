from build_batch6 import faithful
from variants import emit, set_la_thickness_equiv

m = faithful(); set_la_thickness_equiv(m, 1.5); emit('T3f_t15', m)
