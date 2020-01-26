import numpy as np
import scipy

n_units = 5
tau_i = 24 #12
tau_n = 288 #144
c_e = .001
c_ei = .03
a_ee = 14
a_ei = 10
a_en = 4
theta_e = 6
a_ie = 20
a_ii = 8
a_in = .1
theta_i = 5
a_n = 2
beta = 1
p = 2

def frate_func(x, beta):
	return np.sqrt(x/(1 - np.exp(-beta * x)))

def u_der(u, u_tilde, v_tilde, n_tilde, a_ee, a_ei, a_en, theta_e, input_sig, beta):
	return (-u + frate_func(a_ee * u_tilde - a_ei * v_tilde + a_en * n_tilde - theta_e + input_sig, beta))
	# return -u + np.pi * frate_func(a_ee * u_tilde - a_ei * u_tilde + a_en * n_tilde - theta_e + input_sig, beta)

def v_der(v, u, n, a_ie, a_ii, a_in, theta_i, beta):
	return -v + frate_func(a_ie * u - a_ii * v + a_in * n - theta_i, beta)

def n_der(n, a_n, u, p):
	return -n + a_n*(u**p)*(1-n)


def alpha_tilde(all_alphas, unit_i, c_z, N):
	other_alphas_inds = np.arange(N)
	other_alphas_inds = other_alphas_inds[other_alphas_inds != unit_i]

	other_alphas_sum = all_alphas[other_alphas_inds].sum()
	return (all_alphas[unit_i] + c_z * other_alphas_sum)/(1 + c_z * (N - 1))


n_times = 1000
# timep = np.linspace()
all_alpha_time = np.zeros(shape = [3, n_units, n_times])
all_alpha_time[:, :, 0] = np.array([.093, .2, .02])[:, np.newaxis]
all_alpha_tilde_time = np.zeros(shape = [3, n_units, n_times])

unit_as_factor = np.zeros(shape = [n_units, n_times])

# create a sine wave at a certain. freq
freq = 4; start_t = 0; end_t = 2.
x = np.linspace(start_t, end_t, n_times)
sine = np.sin(x*freq*np.pi*2)

for t in np.arange(n_times - 1):
	input_sig = np.zeros(n_units)
	if (t >= 50) & (t < 61):
		input_sig[:1] = 6

	# # add out-of-phase pop
	t_start2 = 200
	if (t >= t_start2) & (t < (t_start2+10)):
		input_sig[2] = 4

	# # # add out-of-phase pop
	t_start2 = 420
	if (t >= t_start2) & (t < (t_start2+10)):
		input_sig[4] = 3

	# # kill pop 2 and 3
	# if (t >= 450) & (t < 600):
	# 	input_sig[4] = 50

	# additional input to pop 2
	# with 3 and 2 S-OP pops you need 500-600 with 20 to shut down the 2-S pops OP with the 3-S pops
	# with 3 and 1 S-OP pops you need 480-580 with 20 to shut down the 1 OP pop
	# if (t >= 500) & (t < 600):
	# 	input_sig[1] = 20

	# compute tildes
	for alpha_ind in np.arange(3):
		if alpha_ind == 1:
			c_z = c_ei
		else:
			c_z = c_e

		for unit_n in np.arange(n_units):
			all_alpha_tilde_time[alpha_ind, unit_n, t + 1] = alpha_tilde(all_alpha_time[alpha_ind, :, t], unit_n, c_z, n_units)

	for alpha_ind in np.arange(3):
		for unit_n in np.arange(n_units):
			u, v, n = all_alpha_time[:, unit_n, t]
			if alpha_ind == 0:
				u = all_alpha_time[alpha_ind, unit_n, t]
				u_tilde, v_tilde, n_tilde = all_alpha_tilde_time[:, unit_n, t]
				all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] +\
					u_der(u, u_tilde, v_tilde, n_tilde, a_ee, a_ei, a_en, theta_e, input_sig[unit_n], beta)

			elif alpha_ind == 1:
				all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] + v_der(v, u, n, a_ie, a_ii, a_in, theta_i, beta)/tau_i

			elif alpha_ind == 2:
				all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] + n_der(n, a_n, u, p)/tau_n

	# send a burst at theta frequency
	if sine[t] >.9:
		for unit_n in np.arange(n_units):
			unit_as_factor[unit_n, t] = all_alpha_time[0, unit_n, t] * np.random.random()*2





cols = ['r', 'g', 'b', 'y', 'k']
fig, axs = pl.subplots(6, 1)
for i in np.arange(5):
	axs[i].plot(all_alpha_time[0, i, :], 'b')
	axs[i].plot(all_alpha_time[1, i, :], 'r')
	axs[i].plot(all_alpha_time[2, i, :], 'k')
	axs[i].set_ylim([-1, 12])
	axs[i].grid()
	axs[-1].plot(unit_as_factor[i, :], color=cols[i])
axs[-1].grid()
axs[-1].plot((sine+1)*10, '--', color='purple')

# for i in np.arange(5):
# 	for t in np.arange(n_times - 1):
# 		print('%i - %.2f' % (t, t % 1000./freq))
# 		# show burst at theta frequency
# 		if (not (t % 200)) & (t > 0):
# 			axs[i].vlines(x = t, ymin=0, ymax=12, linestyles='dashed')



















