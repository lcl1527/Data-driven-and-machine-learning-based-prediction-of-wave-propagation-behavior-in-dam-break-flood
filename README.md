# dam-break-ESN-LSTM
It includes three folders, namely 1D S-V, LSTM and RC-ESN.
(1) 1D S-V: 1D S-V.m is the matlab code used to solve the one-dimensional dam break; 1D S-V data.mat is the solution result; ic0_H.mat is the data set used for training and testing.
(2) LSTM: LSTM.ipynb is the code of LSTM method; lookback4_hd800_nl2_Shift_28.mat is the numerical solution, predicted value, RMSE and ACC of 28 test sets; training set size.xlsx records the predicted horizon and predicted average RMSE under different training set sizes.
(3) RC-ESN: RC-ESN.ipynb is the code of RC-ESN method; Res_size1400_Rd_0.1_Shift_number_28.mat is the numerical solution, predicted value, RMSE and ACC of 28 test sets; Res_size_25_Rd_0.1_Shift_number_1.mat is the RMSE under different reservoir sizes; Res_size_1400_Rd_100_Shift_number_1.mat is the RMSE under different spectral radius; RMSE.xlsx records the predicted horizon and predicted average RMSE under different traininf set sizes.
