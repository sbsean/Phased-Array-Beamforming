% <Orchid BE Simulator Function>
% Scan Conversion
%
% function [sc_data]=scan_converter(input_buffer,register_file_name,...
%                             array_type,interp_mode,scanline_num,sample_number,...
%                             image_x,image_y,...
%                             depth_st,depth_end,convex_radius,angle,c,...
%                             fs,linear_array_length,hi_bound);
%
% sc_data : output data
%
% input_buffer : input data
% sc_coef : sc coeff. array 32x4
% array_type : 0 (Convex Array), 1 (Phased Array), 2 (Linear Array)
% interp_mode : SC interpolation mode : 0(NNI), 1(2x2 linear), 5 (siemens 4x2)
% scanline_num : scanline # for full image
% sample_number : sample #/scanline
% image_x : SC output image width (pixel)
% image_y : SC output image height (pixel)
% depth_st : start depth(mm)
% depth_end : end depth(mm)
% convex_radius : radius of convex array (mm)
% angle : view angle (radian)
% c : ultrasound velocity (mm/s)
% fs : sampling freq.
% linear_array_length : array width (mm) : linear array only
% hi_bound : upper bound for output data saturation

function [sc_data]=scan_converter_rev1(input_buffer,sc_coef,...
                             array_type,interp_mode,scanline_num,sample_number,...
                             image_x,image_y,...
                             depth_st,depth_end,convex_radius,angle,c,...
                             fs,linear_array_length,hi_bound);


%fid=fopen(register_file_name,'r');
%Register_file=fscanf(fid,'%x');
%fclose all;

%sc_coef_file='sc_coef_file_name';
%coef=dump(hex2dec('50000'),hex2dec('5007f'),12,1,2,10,Register_file,sc_coef_file,0);
%filter_coef=reshape(coef,32,4);

filter_coef=sc_coef;

% sc parameters
mode=array_type;
filter_sel=interp_mode;
scanline=scanline_num;
sample_num=sample_number;
image_width=image_x;
image_height=image_y;
view_depth_st=depth_st;
view_depth_end=depth_end;
radius=convex_radius;
view_angle=angle;
ultrasound_speed=c;
sampling_rate=fs;
linear_array_length=linear_array_length;
dz=ultrasound_speed/(2*sampling_rate);
depth=view_depth_end-view_depth_st;
dy=(depth)/image_height;
dx=dy;

in_buffer=input_buffer;
temp=zeros(sample_number,1);
in_buffer=[temp in_buffer temp];

[row col]=size(in_buffer);
in_buffer=reshape(in_buffer,1,row*col);   % comlumn wise

switch mode
 case 0 % convex array
  delta_angle=view_angle/(scanline-1);
  x_offset=(image_width-1)*0.5*dx;
  y_offset=radius*cos(0.5*view_angle);
  theta_offset=(0.5*pi)-(view_angle*0.5);
  r=radius;
 case 1 % phased array
  radius=0;
  r=0;
  delta_angle=view_angle/(scanline-1);
  x_offset=(image_width-1)*0.5*dx;
  y_offset=radius*cos(0.5*view_angle);
  theta_offset=(0.5*pi)-(view_angle*0.5);
 case 2 % linear array
  delta_angle=linear_array_length/(scanline-1);
  x_offset=(dx*image_width-linear_array_length)*0.5;
  y_offset=0.0;
  theta_offset = 0.5*pi;
  r=0;
end

out_buffer=zeros(image_height,image_width);
index_x0=0:1:image_width-1;
index_y0=0:1:image_height-1;

% 2d array (x-index, y-index)
index_x=ones(image_height,1)*index_x0;		% row wise
index_y=index_y0'*ones(1,image_width);		% row wise

%[index_x,index_xx]=finargsz('all',index_x0,out_buffer);
%[index_y,index_yy]=finargsz('all',index_y0',out_buffer);

dis_x=index_x*dx-x_offset;
dis_y=index_y*dy+y_offset;

switch mode
 case 0
  s_temp=((0.5*pi)+(atan(dis_x./dis_y))-theta_offset).*(dis_y~=0);
  s=s_temp;
  s_temp=((0.5*pi)+(0.5*pi)-theta_offset).*(dis_y==0);
  s=s+s_temp;
      
  d=sqrt(dis_x.*dis_x+dis_y.*dis_y)-r;
  s=s/delta_angle;
  d=d/dz;
  
 case 1
  s_temp=((0.5*pi)+(atan(dis_x./dis_y))-theta_offset).*(dis_y~=0);
  s=s_temp;
  s_temp=((0.5*pi)+(0.5*pi)-theta_offset).*(dis_y==0);
  s=s+s_temp;
      
  d=sqrt(dis_x.*dis_x+dis_y.*dis_y)-r;
  s=s/delta_angle;
  d=d/dz;
  
 case 2
  s=dis_x;
  d=dis_y;
  s=s/delta_angle;
  d=d/dz;	
end	   

if (filter_sel==0)
	sn=round(s);
	dn=round(d);
else         
	sn=floor(s);
	dn=floor(d);
end


sn=reshape(sn,1,image_height*image_width);
dn=reshape(dn,1,image_height*image_width);
d=reshape(d,1,image_height*image_width);
s=reshape(s,1,image_height*image_width);

data1=zeros(1,image_height*image_width);
data2=zeros(1,image_height*image_width);
data3=zeros(1,image_height*image_width);

% origin : sn+1, dn+1
switch filter_sel
 case 0
  flag=ones(1,image_height*image_width);
  flag=(flag).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
   
  index1=(dn+1)+((sn+2-1)*sample_number);
   
  index1=(index1-1).*(index1<=length(in_buffer))+1;
  index1=(index1-1).*(index1>0)+1;
    
  data1=in_buffer(index1);
  data3=flag.*data1;

 case 1
  flag=ones(1,image_height*image_width);
  flag=(flag).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
% sn=sn.*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
% s=s.*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
% dn=dn.*(dn>=0 & dn<sample_num-1 & sn>=0 & sn<scanline-1);
  
  index1=(dn+1)+((sn+2-1)*sample_number);
  index2=(dn+1)+((sn+3-1)*sample_number);
  index3=(dn+2)+((sn+2-1)*sample_number);
  index4=(dn+2)+((sn+3-1)*sample_number);
  
  index1=(index1-1).*(index1<=length(in_buffer))+1;
  index2=(index2-1).*(index2<=length(in_buffer))+1;
  index3=(index3-1).*(index3<=length(in_buffer))+1;
  index4=(index4-1).*(index4<=length(in_buffer))+1;
  
  index1=(index1-1).*(index1>0)+1;
  index2=(index2-1).*(index2>0)+1;
  index3=(index3-1).*(index3>0)+1;
  index4=(index4-1).*(index4>0)+1;
  
  data1=(sn+1-s).*in_buffer(index1)+(s-sn).*in_buffer(index2);
  data2=(sn+1-s).*in_buffer(index3)+(s-sn).*in_buffer(index4);
  data3=flag.*((dn+1-d).*data1+(d-dn).*data2);
 case 2
 case 3
 case 4
 case 5
%  phase=((floor((s-sn)*32.0)+1)-1).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1)+1;					
%  q=floor((d-dn)*16).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
%  range=q/16.0;
%  flag=ones(1,image_height*image_width);
%  flag=(flag).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
%  sn=(sn+0).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1)-0;
%  s=(s).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
%  dn=dn.*(dn>=0 & dn<sample_num-1 & sn>=0 & sn<scanline-1);

  phase=((floor((s-sn)*32.0)+1)-1).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1)+1;					
  q=floor((d-dn)*16).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
  range=q/16.0;
  flag=ones(1,image_height*image_width);
  flag=(flag).*(sn>=0 & sn<scanline-1 & dn>=0 & dn<sample_num-1);
    
  index1=(dn+1)+((sn+1-1)*sample_number);
  index2=(dn+1)+((sn+2-1)*sample_number);
  index3=(dn+1)+((sn+3-1)*sample_number);
  index4=(dn+1)+((sn+4-1)*sample_number);
  index5=(dn+2)+((sn+1-1)*sample_number);
  index6=(dn+2)+((sn+2-1)*sample_number);
  index7=(dn+2)+((sn+3-1)*sample_number);
  index8=(dn+2)+((sn+4-1)*sample_number);
  
  index1=(index1-1).*(index1<=length(in_buffer))+1;
  index2=(index2-1).*(index2<=length(in_buffer))+1;
  index3=(index3-1).*(index3<=length(in_buffer))+1;
  index4=(index4-1).*(index4<=length(in_buffer))+1;
  index5=(index5-1).*(index5<=length(in_buffer))+1;
  index6=(index6-1).*(index6<=length(in_buffer))+1;
  index7=(index7-1).*(index7<=length(in_buffer))+1;
  index8=(index8-1).*(index8<=length(in_buffer))+1;
  
  index1=(index1-1).*(index1>0)+1;
  index2=(index2-1).*(index2>0)+1;
  index3=(index3-1).*(index3>0)+1;
  index4=(index4-1).*(index4>0)+1;
  index5=(index5-1).*(index5>0)+1;
  index6=(index6-1).*(index6>0)+1;
  index7=(index7-1).*(index7>0)+1;
  index8=(index8-1).*(index8>0)+1;
  
  whos filter_coef(phase,4)
  whos in_buffer(index1)
  
  data1=filter_coef(phase,4)'.*in_buffer(index1)+filter_coef(phase,3)'.*in_buffer(index2)+...
        filter_coef(phase,2)'.*in_buffer(index3)+filter_coef(phase,1)'.*in_buffer(index4);

  data2=filter_coef(phase,4)'.*in_buffer(index5)+filter_coef(phase,3)'.*in_buffer(index6)+...
        filter_coef(phase,2)'.*in_buffer(index7)+filter_coef(phase,1)'.*in_buffer(index8);
  
  data3=flag.*((1-range).*data1+(range).*data2);

end

data3=reshape(data3,image_height,image_width);
out_buffer=floor(data3);
out_buffer=out_buffer.*(out_buffer>=0)+0;
out_buffer=(out_buffer-hi_bound).*(out_buffer<=hi_bound)+hi_bound;
sc_data=out_buffer;

%im=uint8(out_buffer);
%imwrite(im,out_file,'bmp');




	