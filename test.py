from bin_option import *

option_values = [50,50,0.04,0.4,(183/365),0.01]
simulations = 50000
bs_text = ['(0.46473592872321867, 0.53526407127678133)','(28.715375330447621, 21.03456711845768)']
OPTION_ITEMS = 'Black Scholes', 'Binomial', 'Monte-Carlo', 'Implicit', 'Explicit', 'Crank-Nicolson'

def run_test_output():
    print('Testing for: ' + bo.bin_option_types[bo.bin_asset_type-1])
    print('Black Scholes')
    print('E: ' + bs_text[bo.bin_asset_type-1])
    print('A: ' + str(bo.binary_black_scholes()))
    print('Monte-Carlo')
    print(bo.binary_monte_carlo(simulations))

    print('=====')

def set_type(new_type):
    new_type = str(new_type)
    print('Setting: ' + new_type)
    try:
        bo.set_bin_asset_type(new_type)
        print('Okay')
    except:
        print('Cannot')
    print('=====')

def get_z(N):
    z = np.random.standard_normal(N)
    z = np.append(z,-z)
    return z

bo = BinaryOption()
bo.set_variables(option_values)
outputs = {
        0: bo.black_scholes(),
        1: bo.binomial_tree(),
        2: bo.monte_carlo(get_z(10_000_000)),
        3: bo.implicit(),
        4: bo.explicit(),
        5: bo.crank_nicolson(),
    }

run_test_output()
set_type(2)
run_test_output()
set_type('cash')
set_type('asset')
set_type('vanilla')
set_type(5)
bo.set_bin_asset_type(2)
print('Finally: ' + bo.bin_option_types[bo.bin_asset_type-1])
print('=====')
for i,v in enumerate(outputs):
    print('Testing ' + OPTION_ITEMS[i])
    print(outputs[i])
    print('=====')