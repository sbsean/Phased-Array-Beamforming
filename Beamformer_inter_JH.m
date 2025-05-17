clc; clear all; close all;
format long g

Beam_Num=96;				            % scanline #
Elem_Pitch=0.26;		            	% element pitch [mm]
Elem_Num=64;			            	% element #
Ch_Num=16;				                % active Ch. #
Sound_Speed=1540000;	            	% [m/s]
Master_Clock=40e6;			            % [Hz]
Unit_Dis=Sound_Speed/Master_Clock;	    % [mm]
Delta_z=Sound_Speed/(2*Master_Clock);	% [mm]
View_Angle=80;                          % 80 degree view 
Delta_Angle=View_Angle/(Beam_Num-1);	% [degree]
data_length=10240;

%파일로드 형식 RxScanline0.bin, RxScanline1.bin, RxScanline2.bin, ... ,
%RxScanline95.bin
%Beamformed_data 폴더에 sum_out.bin으로 저장