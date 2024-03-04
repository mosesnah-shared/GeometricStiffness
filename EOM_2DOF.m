% ==============================================================================
%   Script Name: EOM_2DOF.m
%        Author: Moses C. Nah
%         Email: mosesnah@mit.edu
%   Institution: MIT
%  Date Created: 03-04-2024
% Last Modified: 03-04-2024
% 
% Description:
%   This script derives the equations of motion of a 2DOF robot
%   with revolute joints.
%
% ==============================================================================

%% (--) Initialization
clear; close all; clc;

%% (1A) Euler-Lagrange Equations

% Geometric and Inertial Parameters 
syms M1 M2 L1 L2 positive

% Relative Angle Coordinates
syms t 
syms q1(t) q2(t)
q_arr = { q1, q2 };

% The Kinematics of both point masses
x1 = L1*cos( q1 );
y1 = L1*sin( q1 );

x2 = L1*cos( q1 ) + L2*cos( q1+q2 );
y2 = L1*sin( q1 ) + L2*sin( q1+q2 );

dx1 = diff( x1, t ); dy1 = diff( y1, t );
dx2 = diff( x2, t ); dy2 = diff( y2, t );
dq1 = diff( q1, t ); dq2 = diff( q2, t );

% Lagrangian == Total Kinetic Energy of the System
L = 1/2*M1*( dx1^2 + dy1^2 ) + 1/2*M2*( dx2^2 + dy2^2 );

% The Equations of motion
eq1 = simplify( diff( diff( L, dq1 ), t ) - diff( L, q1 ) );
eq2 = simplify( diff( diff( L, dq2 ), t ) - diff( L, q2 ) );

% Getting the Mass matrix 
m11 =  diff( diff( L, dq1 ), dq1 );
m12 =  diff( diff( L, dq1 ), dq2 );
m22 =  diff( diff( L, dq2 ), dq2 );

M    = simplify( [ m11, m12; m12, m22 ] );
Minv = simplify( inv( M ) );

% Substitute the Symbolic Form
M    = subs(    M, { q1, q2 } , { 'q1', 'q2' } );
Minv = subs( Minv, { q1, q2 } , { 'q1', 'q2' } );


%% (1B) The Levi-Civita Connections

Gamma = sym( zeros( 2, 2, 2 ) );

for i = 1:2
    for j = 1:2
        for k = 1:2
           tmp = 0;
           for m = 1:2
               tmp = tmp + 1/2*Minv( k, m ) * ( diff( M( m, i ), q_arr{ j } ) + diff( M( m, j ), q_arr{ i } ) - diff( M( i, j ), q_arr{ m } ) );
           end
           Gamma( i, j, k ) = tmp;
        end
    end
end