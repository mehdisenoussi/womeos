from psychopy import visual, event, core
import os, time, glob
import numpy as np
import pylab as pl
from psychopy.visual import ShapeStim

# Create the Psychopy window
myWin = visual.Window(fullscr=False, screen=0, monitor='testMonitor', allowGUI=True) #, mouseVisible=False)
myWin.setMouseVisible(True)


frame_period = myWin.monitorFramePeriod
send_eegtriggers = False

####################################################
# Create Psychopy visual objects
####################################################
fixcross_size=.15
fixh = visual.Line(myWin, units='deg', start=(-fixcross_size, 0), end=(fixcross_size, 0), lineColor=(-1,-1,-1), lineWidth=2)
fixv = visual.Line(myWin, units='deg', start=(0, -fixcross_size), end=(0, fixcross_size), lineColor=(-1,-1,-1), lineWidth=2)

# create shapes textures
shapes_sizes = np.array([2, 2, 2, 1])
positions = np.array([[-2, -2], [-2, 2], [2, -2], [2, 2]])
shapes_textures = np.array([	visual.ShapeStim(myWin, units='deg', lineWidth=6, lineColor='black', vertices=((0, .6), (.4, -.2), (-.4, -.2)),
									pos=positions[0], fillColor='white', size=shapes_sizes[0]),#									  1		     2		   3			4		  5		 	 6		     7		 		8		 9			 10
								visual.ShapeStim(myWin, units='deg', lineWidth=6, lineColor='black', vertices=((0, .6), (.15, .25), (.6, .25), (.25, 0), (.41, -.45), (0, -.18), (-.41, -.45), (-.25, 0), (-.6, .25), (-.15, .25)),
									pos=positions[2], fillColor='white', size=shapes_sizes[1]),
								visual.ShapeStim(myWin, units='deg', lineWidth=6, lineColor='black', vertices=((-.4, .4), (.4, .4), (.4, -.4), (-.4, -.4)),
									pos=positions[1], fillColor='white', size=shapes_sizes[2]),
								visual.Polygon(myWin, edges=5, units='deg', lineWidth=6, lineColor='black', radius=shapes_sizes[3], pos=positions[3], fillColor='white')])


timer = core.CountdownTimer()

# colors: blue, green, red, yellow
colors = np.array([[38, 160, 218], [104, 227, 102], [249, 74, 43], [255, 223, 43]])
col_names = np.array(['blue', 'green', 'red', 'yellow'])

# create cue text texture
lett_size = 1; lett_pos = 1
cue_text_textures = [visual.TextStim(myWin, units='deg', height = lett_size,
									pos=(0, lett_pos), text=col, alignHoriz = 'center', color='black') for col in col_names]

n_trials = 3

loads = np.arange(1, 6)
col_inds = np.arange(4)
sizes = np.arange(3)
timings = np.zeros(shape = [200, 7])

# minimum number of levels across all manipulable dimensions
min_n_levels = np.min([len(loads), len(col_inds), len(sizes)])
# how many times do we have to repeat the minimal number of levels to get to n_trials trials
n_reps_min = int(np.ceil(n_trials/min_n_levels))

# make an array containing all trial info
trial_info = np.array([np.repeat(loads, n_reps_min)[:n_trials], np.tile(col_inds, n_reps_min)[:n_trials], np.repeat(sizes, n_reps_min)[:n_trials]])

prestimdur = 1.
stim_dur = .5
delay_1 = 1.
cue_dur = .5
delay_2 = 1.
target_dur = 1.
ISIdur = 1.


txt=''; trial_n=0
while txt!='Bye!' and trial_n < n_trials:
	####### trial baseline (pre-instruction) #######
	fixh.lineColor=[-1,-1,-1]; fixv.lineColor=[-1,-1,-1]
	fixh.draw(); fixv.draw()
	timings[trial_n,0] = myWin.flip()
	timer.reset(prestimdur - (frame_period/2.))
	if send_eegtriggers:
		s.write('A010')
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()

	####### Stim array #######
	# parametrize and draw stimuli
	rand_colors = np.random.choice(len(colors), 4, replace=True)
	rand_size = np.random.choice([.9, 1, 1.1], 4, replace=True)
	rand_pos = np.random.choice(len(shapes_textures), 4, replace=False)
	
	for shape_i in np.arange(len(shapes_textures)):
		shapes_textures[shape_i].fillRGB = colors[rand_colors[shape_i], :]
		shapes_textures[shape_i].size = np.array([rand_size[shape_i], rand_size[shape_i]]) * shapes_sizes[shape_i]
		shapes_textures[shape_i].pos = positions[rand_pos[shape_i], :]
		shapes_textures[shape_i].draw()
	
	fixh.draw(); fixv.draw()

	while timer.getTime()>0:
		continue
	timings[trial_n,1] = myWin.flip()

	timer.reset(stim_dur - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(20+trials_info[trial_n,0]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()
	
	# while timer.getTime() > 0:
	# 	continue

	# trial_n+=1

	####### Delay 1 #######
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 2] = myWin.flip()
	timer.reset(delay_1 - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()

	####### Cue #######
	cued_color_ind = np.random.randint(len(col_names))
	cue_text_textures[cued_color_ind].draw()
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 3] = myWin.flip()
	timer.reset(cue_dur - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()

	####### Delay 2 #######
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 4] = myWin.flip()
	timer.reset(delay_2 - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()


	####### Target #######
	rand_ori = np.random.choice([-10, 0, 10], 1, replace=True)[0]
	rand_target = np.random.randint(4)
	rand_size = np.random.choice([.9, 1, 1.1], 1, replace=True)
	shapes_textures[rand_target].pos = np.array([0, 0])
	shapes_textures[rand_target].size = (np.array([rand_size, rand_size]) * shapes_sizes[rand_target]).squeeze()
	shapes_textures[rand_target].ori = rand_ori
	shapes_textures[rand_target].fillRGB = None
	shapes_textures[rand_target].draw()
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 5] = myWin.flip()
	timer.reset(delay_1 - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()


	####### ISI #######
	fixh.lineColor=[1,1,1]; fixv.lineColor=[1,1,1]
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 6] = myWin.flip()
	timer.reset(ISIdur - (frame_period/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()
	while timer.getTime() > 0:
		continue

	trial_n += 1












if train_settings:
	instdur = 1
	ISIdur =  np.arange(3,3.5,.05)
	stimdur = 1
	max_respTime = 1.5


# Create trial structure (instructions, ISIs and stim combination (Hori-Vert, Hori-Hori, etc.))
trials_info = exputils.make_mystinfo_design(n_instructions=4, n_ISIs=ISIdur.shape[0], n_stimCombi=4, block=block, hash_ISIs_by_block=True)
trials_info = trials_info.astype(np.uint8)
n_trials = trials_info.shape[0]
# Create key bindings
if finger_asso == 'a':
	key_asso = {'vert_keys':['s', 'k'], 'hori_keys':['d', 'l'], 'left_keys':['s', 'd'], 'right_keys':['k', 'l']}
elif finger_asso == 'b':
	key_asso = {'vert_keys':['s', 'l'], 'hori_keys':['d', 'k'], 'left_keys':['s', 'd'], 'right_keys':['k', 'l']}
resp_keys = np.hstack([key_asso['vert_keys'], key_asso['hori_keys']])
all_possible_keys = np.hstack([key_asso['vert_keys'], key_asso['hori_keys'], 'escape', 'esc', 'space'])


####################################################
# Show instructions for the experiment
####################################################
# clear the back buffer for drawing
myWin.clearBuffer()

slide_n = 0
while slide_n < 7:
	if (slide_n + 1) in [2, 4, 5]: suffix = finger_asso
	else: suffix = ''
	instr_protocol=visual.ImageStim(myWin, './material/protocol_instruction_slides/slide_%i' % (slide_n+1) + suffix + '.png')
	instr_protocol.draw(); myWin.flip()
	k = event.waitKeys(keyList=['right', 'left', 'space'])

	if k[0] == 'left': slide_n -= 1
	elif k[0] == 'right' or k[0] == 'space': slide_n += 1

	if slide_n < 0: slide_n = 0

myWin.flip()
time.sleep(3)
####################################################



# for trial_n in range(1):
txt=''; trial_n=0
while txt!='Bye!' and trial_n < n_trials:

	####### trial baseline (pre-instruction) #######
	fixh.lineColor=[-1,-1,-1]; fixv.lineColor=[-1,-1,-1]
	fixh.draw(); fixv.draw()
	timings[trial_n,0] = myWin.flip()
	timer.reset(preinstdur - (ref_rate/2.))
	if send_eegtriggers:
		s.write('A010')
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()

	####### Instruction #######
	# choose which letters to display for hand to use and target stimulus
	# order of instructions (hand-stim): 0 = L-L, 1 = L-R, 2 = R-L, 3 = R-R

	# hand letter is L if instruction index is 0 or 1, R if instruction index is 2 or 3
	hand_tex = hand_lett_textures[int(trials_info[trial_n,0]>1)]
	# stim letter is L if instruction index is even (0 or 2), R if odd (1 or 3)
	stim_tex = stim_lett_textures[int(trials_info[trial_n,0]%2)]

	# draw instructions
	hand_tex.draw(); stim_tex.draw()
	fixh.draw(); fixv.draw()

	while timer.getTime()>0:
		continue
	timings[trial_n,1] = myWin.flip()

	timer.reset(instdur - (ref_rate/2.))
	if send_eegtriggers:
		s.write('A0%i'%(20+trials_info[trial_n,0]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()



	####### ISI #######
	fixh.draw(); fixv.draw()
	ISI_trial_n[trial_n] = ISIdur[trials_info[trial_n, 1]]
	while timer.getTime() > 0:
		continue

	timings[trial_n,2] = myWin.flip()
	timer.reset(ISI_trial_n[trial_n] - (ref_rate/2.))
	if send_eegtriggers:
		s.write('A0%i'%(30 + trials_info[trial_n, 1]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()


	####### Grating #######
	if trials_info[trial_n,2] == 0:
		myStim[0].ori, myStim[1].ori=[0, 0]
	elif trials_info[trial_n,2] == 1:
		myStim[0].ori, myStim[1].ori=[0, 90]
	elif trials_info[trial_n,2] == 2:
		myStim[0].ori, myStim[1].ori=[90, 0]
	elif trials_info[trial_n,2] == 3:
		myStim[0].ori, myStim[1].ori=[90, 90]

	myStim[0].draw(); myStim[1].draw()
	fixh.draw(); fixv.draw()

	while timer.getTime() > 0:
		continue
	timings[trial_n, 3] = myWin.flip()
	temp_resp_t = time.time()
	timer.reset(stimdur - (ref_rate/2.));
	if send_eegtriggers:
		s.write('A0%i'%(50 + trials_info[trial_n, 2]))
		time.sleep(.004)
		s.write('A000')
	myWin.clearBuffer()

	fixh.lineColor=[-1, -1, 1]; fixv.lineColor=[-1, -1, 1];
	fixh.draw(); fixv.draw()
	while timer.getTime() > 0:
		continue
	timings[trial_n, 4] = myWin.flip()

	k = event.waitKeys(maxWait=max_respTime, keyList=all_possible_keys, timeStamped=timings[trial_n, 3])
	resp_times2 = time.time() - temp_resp_t

	resp_key='noresp'
	correct=False
	if k is None:
		txt = 'you did NOT respond'
		k =[['none', -1.]]
		wtxt = 'Too late! Do you need a break? press SPACE to continue..'
		waitText = visual.TextStim(myWin, units='pix',height = 20,
			pos=(-400, 0), text=wtxt, alignHoriz = 'left', color='black')
		waitText.draw(); myWin.flip()
		kk = event.waitKeys(keyList=['space'])

		# add the trial to the end of the list
		trials_info=np.concatenate([trials_info, trials_info[trial_n,:][np.newaxis]], axis=0)
		n_trials+=1
		resp_times[trial_n] = -1;

	elif k[0][0] in resp_keys:
		# if observer responded, check correctness and change fixation cross color
		resp_key=k[0][0]
		resp_times[trial_n] = k[0][1]
		correct = exputils.compute_resp_correctness(key=resp_key, trial_info=trials_info[trial_n,:], key_asso=key_asso)
		corr_resps[trial_n] = correct
		if send_eegtriggers:
			s.write('A0%i'%(60+correct))
			time.sleep(.004)
			s.write('A000')
			
		rt_txt = '%i' % int(resp_times[trial_n]*1000)
		resptimeText = visual.TextStim(myWin, units='pix',height = 20,
			pos=(0, -20), text=rt_txt, alignHoriz = 'center', color='black')
		if correct: fixh.lineColor=[-1,1,-1]; fixv.lineColor=[-1,1,-1];
		else: fixh.lineColor=[1,-1,-1]; fixv.lineColor=[1,-1,-1];
		fixh.draw(); fixv.draw()
		resptimeText.draw()
		myWin.flip()
		timer.reset(.5 - (ref_rate/2.))
		myWin.clearBuffer()
		fixh.lineColor=[1,1,1]; fixv.lineColor=[1,1,1];
		fixh.draw(); fixv.draw()
		while timer.getTime()>0:
			continue

	elif k[0][0] in ['escape', 'esc']: break

	task_out.write('%s\t%i\t%i\t%i\t%i\t%i\t%s\t%.5f\t%.5f\t%i\n'%(obs,\
		block, trial_n, trials_info[trial_n,0], trials_info[trial_n,1], trials_info[trial_n,2], resp_key, resp_times[trial_n], resp_times2, correct))


	####### ITI #######
	timings[trial_n, 5] = myWin.flip()
	timer.reset(np.random.choice(ITIdur) - (ref_rate/2.))
	myWin.clearBuffer()
	while timer.getTime()>0:
		continue

	trial_n+=1

# send end of block trigger
if send_eegtriggers: s.write('A100'); time.sleep(.004); s.write('A000')

# close the log file
task_out.close()

# save parameters
outfolder = './results/obs_%s/' % obs
output_param_file = '%sparameters_obs%s_block%i_date%s.npz' % (outfolder, obs, block, timeAndDate)
np.savez(output_param_file, {'trials_info':trials_info, 'timings':timings, 'resp_times':resp_times, 'stimSize':stimSize, 'contrast':contrast,\
		'oris':oris, 'dist_to_center':dist_to_center, 'spatFreq':spatFreq, 'ref_rate':ref_rate, 'preinstdur':preinstdur,\
		'instdur':instdur, 'ISIdur':ISIdur, 'stimdur':stimdur, 'ITIdur':ITIdur, 'max_respTime':max_respTime} )


if k[0][0]=='escape': txt='ESCAPED ! bye !'
else: txt='FINISHED ! Thanks !'
endText = visual.TextStim(myWin, units='pix', height = 30,
            pos=(0, 0), text=txt, alignHoriz = 'center', color='black')
endText.draw(); myWin.flip()
time.sleep(2)

print('\n\n\tAccuracy on block %i was: %.1f\n\n' % (block, corr_resps[:n_trials].mean()))

myWin.close()
core.quit()


