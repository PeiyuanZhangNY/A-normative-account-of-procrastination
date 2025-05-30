# A-normative-account-of-procrastination
# Files in the folder Preprint figure:
Program for Figure Model and functions: rewardcostfunction.m; Final figure: combined_modelfigure_updated.pdf
Program for Figure Temporal patterns of work progress and its associated parameter space: pattern3D.m; Final figure: pattern3Drevised.pdf
Program for Figure Effects of maximum task reward, task aversion, and utility of alternative activity: alphatimecoursecurve_daysofpro_utility.m; Cmaxtimecoursecurve_daysofpro_utility.m; Jtimecoursecurve_daysofpro_utility.m; Final figure: combined_alpha_Cmax_J_allutility.pdf
Program for Figure Effects of discount rates: gammatimecoursecurve_daysofpro_utility.m; Final figure: gamma_allutility_plot.pdf
Program for Figure Effects of the shape of the cost function: lambdatimecoursecurve_daysofpro_utility.m; Final figure: lambda_allutility_plot.pdf
Program for Figure Diverse relationship between perfectionism and procrastination: Betatimecoursecurve_daysofpro_utility.m; Final figure: Betatimecourse_allgamma_allutility.pdf
Program for Figure Effects of total given time: Ttimecoursecurve_daysofpro_utility.m; Final figure: T_allutility_plot.pdf
Draft for Figure Illustrations of interventions (Please request access to it. Correspondence to Pei Yuan Zhang, Email: pz580@nyu.edu): https://docs.google.com/presentation/d/1nN4_nmuHwHJH2gKKn416oGMjkEy8Er8j6NHz0JLBAm8/edit#slide=id.g275b3ee0b89_0_0
Program for Figure Effects of all immediate reward interventions: 
PatternTimeCourse_DelayedReward_ImmediateRewardInterventions.m 
combinedPatternTimeCourseRewardMilestoneEachUnitProgress.m
daysofpro_DelayedReward_ImmediateRewardInterventions.m; 
Final figure: immediate_reward_interventions_withdiagram_4.pdf
Program for Figure Effects of intermediate deadline interventions: interim_deadline_intervention.m; Final figure: interim_deadline_interventions_withdiagram revised_2.pdf

# Files in the folder subfigures: 
Program for Figure Supplementary Figure 2: perfectionismswitchNotWorkingtoCompleting.m; Final figure: perfectionismswitchNotWorkingtoCompleting.pdf

# Files in the folder Supplement:
DelayedRewardGammaPatternProcrastination.m simulation program supporting the necessary parameter space where we can find “a delay in the beginning, and then ramping up. ” The results saved in DelayedRewardGammaPatternProcrastination.mat

DelayedRewardConcaveCost.m simulation program supporting the parameter space where we can find “working at the last minute”. The results saved in DelayedRewardConcaveCost.mat

necessaryconditionNotWorking.m simulation program supporting the necessary parameter space where we can find “not working at all”. The results saved in necessaryconditionNotWorking.mat 

alphaMonteCarlo.m simulation program supporting the finding: under a lower maximum reward, agents procrastinate for longer, finish less work in the end, and pay a lower total cost. The results are saved in alphaMonteCarlo.mat

CmaxMonteCarlo.m simulation program supporting the finding about how the time course of progress, days of procrastination, final proportion completed, and total cost change under various maximum cost. The results for the non-zero utility of alternative activities are saved in CmaxMonteCarlo_Jnonzero.mat and the results for zero utility of alternative activities are saved in CmaxMonteCarlo_Jzero.mat. 

JMonteCarlo.m simulation program supporting the finding about how the time course of progress, days of procrastination, final proportion completed, and total cost change under various utility of alternative activity. The results are saved in JMonteCarlo.mat

GammaMonteCarlo.m simulation program supporting the finding about how the time course of progress, days of procrastination, final proportion completed, and total cost change under various discount rates. The results are saved in GammaMonteCarlo.mat

lambdaMonteCarlo.m simulation program supporting the finding about how the time course of progress, days of procrastination, final proportion completed, and total cost change under various exponent of cost function. The results for non-zero utility of alternative activities are saved in lambdaMonteCarlo_Jnonzero.mat and for zero utility of alternative activities are saved in lambdaMonteCarlo_Jzero.mat. 

perfectionismFinalCompletion.m simulation program supporting the following finding: For agents with high perfectionism, the final proportion completed is either 0 or 1. The results are saved in perfectionismFinalCom_betalargerlambda.mat. For agents with low perfectionism, the final proportion completed is between 0 and 1, including 0 and 1. The results are saved in perfectionismFinalCom_betasmallerlambda.mat

TMonteCarlo.m simulation program supporting the following finding: Given a longer total time, agents switch from not working at all to working all day and finish more work in the end (we only discuss the case of a convex cost function). The change in the total cost is bidirectional. The results are saved in TMonteCarlo.mat

The following programs are simulations for hyperbolic discounting:  
•	Hyperbolic_GammaPatternProcr.m simulation program supporting the necessary parameter space where we can find “a delay in the beginning, and then ramping up. ” It uses the function program HyperbolicDiscountingFun.m and Qfun.m. The results are saved in Hyperbolic_GammaPatternProcr.mat.
•	Hyperbolic_ConcaveCost.m simulation program supporting the parameter space where we can find “working at the last minute”. The results are saved in Hyperbolic_ConcaveCost.mat
•	Hyperbolic_necessaryconditionNotWorking.m simulation program supporting the necessary parameter space where we can find “not working at all”. The results are saved in Hyperbolic_necessaryconditionNotWorking_Jsmaller.mat and Hyperbolic_necessaryconditionNotWorking_Jlarger.mat
•	Hyperbolic_alphaMonteCarlo.m simulation program supporting the finding: under a lower maximum reward, agents procrastinate for longer, finish less work in the end. The results are saved in Hyperbolic_alphaMonteCarlo.mat. 
•	Hyperbolic_CmaxMonteCarlo.m simulation program for maximum cost. The results are saved in Hyperbolic_CmaxMonteCarlo_Jzero.mat and Hyperbolic_CmaxMonteCarlo_Jnonzero.mat.
•	Hyperbolic_JMonteCarlo.m simulation program for the utility of alternative activity. The results are saved in Hyperbolic_JMonteCarlo.mat.
•	Hyperbolic_k_DRMonteCarlo.m simulation program for hyperbolic discount rate. The results are saved in Hyperbolic_k_DRMonteCarlo.mat.
•	Hyperbolic_perfectionismFinalCompletion.m simulation program for findings of perfectionism. The results are saved in hyperbolic_perfectionism_betalargerlambda.mat and hyperbolic_perfectionism_betasmallerlambda.mat.
•	Hyperbolic_TMonteCarlo.m simulation program for findings of number of total given days. The results are saved in Hyperbolic_TMonteCarlo.mat. 

The following programs are intervention simulations for hyperbolic discounting: 
•	immediateRewardInterventions_hyperbolicdiscount.m for immediate reward interventions. The results are saved in immediateRewardInterventions_hyperbolicdiscount.mat. 
•	InterimDeadlineInterventions_hyperbolicdiscount.m for interim deadline interventions. The results are saved in InterimDeadlineInterventions_hyperbolicdiscount.mat. 


# Files in the folder perfectionism data and analysis:
perfectionism_rawdata.xlsx : The raw data of participants’ response to Self-oriented subscale of perfectionism from Hewitt, P. L., & Flett, G. L. (1990) and response to our designed questionnaire to test a power-law relationship between reward and final performance. 

perfectionism_model.m Model fitting program for a power-law relationship between reward and final performance.
perfectionism_DataAnalysis.ipynb data analysis program for the section: Empirical evidence supporting our model’s potential to represent perfectionism. 
Final figures: satisfactionlevel_modelfit.pdf; perfectionism correlation.pdf

