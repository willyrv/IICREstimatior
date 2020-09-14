def stationary_nisland_iicr(m, n):
    from math import sqrt, exp

    alpha = m * n + n - 1
    delta = sqrt(alpha * alpha - 4 * m * (n - 1))
    coef_e1 = (delta - alpha) / (2 * (n - 1))
    coef_e2 = (-delta - alpha) / (2 * (n - 1))
    coef_f1 = (delta - alpha + 2 * m) / (2 * delta)
    coef_f2 = m * (n - 1) / delta 

    def iicr(t):
        e1 = exp(t * coef_e1)
        e2 = exp(t * coef_e2)
        f1 = coef_f1 * e1 + (1 - coef_f1) * e2
        f2 = coef_f2 * (e1 - e2)
        
        d_e1 = coef_e1 * e1
        d_e2 = coef_e2 * e2
        d_f1 = coef_f1 * d_e1 + (1 - coef_f1) * d_e2
        d_f2 = coef_f2 * (d_e1 - d_e2)

        return -(f1 + f2) / (d_f1 + d_f2)
    
    return iicr

def generate_panmictic_ms_from_data(time, iicr):
	import subprocess
	import numpy as np
	
	ms_seq_simulations = 100
	ms_theta = 1
	ms_recombination = 3
	ms_sites = int(1e4)

	ms_command = ['./scrm', '2', str(ms_seq_simulations), 
		'-t', str(ms_theta), '-r', str(ms_recombination),
		str(ms_sites), '-p', '8']
	
	for i in range(len(time)):
		ms_command.extend(['-eN', str(0.5 * time[i]), str(iicr[i])])

	# run scrm command
	msout_filename = 'panmictic_ms_output.scrm'
	msout = open(msout_filename, 'w')
	proc = subprocess.Popen(ms_command, stdout = msout)
	proc.wait()
	msout.close()

def generate_panmictic_ms_from_parameters(migration_rate, islands, time_limits = [-2, 2], samples = 64):
	import numpy as np
	
	iicr_function = stationary_nisland_iicr(migration_rate, islands)

	time = [10**t for t in np.linspace(*time_limits, samples)]
	iicr = [iicr_function(t) for t in time]
	generate_panmictic_ms_from_data(time, iicr)
