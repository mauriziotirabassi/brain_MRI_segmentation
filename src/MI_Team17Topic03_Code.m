%% DATA EXTRACTION
clear variables; close all; clc

% Setting the path to the resources relatively to the running OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-mi-2023/src';
else
    path = '';
end