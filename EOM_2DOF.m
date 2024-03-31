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
eq1 = simplify( diff( diff( L, dq1 ), t ) - diff( L, q1 ) ); eq1 = eq1( t );
eq2 = simplify( diff( diff( L, dq2 ), t ) - diff( L, q2 ) ); eq2 = eq2( t );

% Getting the Mass matrix 
m11 =  diff( diff( L, dq1 ), dq1 );
m12 =  diff( diff( L, dq1 ), dq2 );
m22 =  diff( diff( L, dq2 ), dq2 );

M    = simplify( [ m11, m12; m12, m22 ] );
Minv = simplify( inv( M ) );

% Taking off the (t) variable in q1, q2 and also from M
M    = subs(    M, { q1, q2 }, { 'q1', 'q2' } );    M = M( t );
Minv = subs( Minv, { q1, q2 }, { 'q1', 'q2' } ); Minv = Minv( t );

% The Levi-Civita Connections
Gamma = sym( zeros( 2, 2, 2 ) );
q_arr = { 'q1', 'q2' };

for i = 1:2
    for j = 1:2
        for k = 1:2
           tmp = 0;
           for m = 1:2
               tmp = tmp + 1/2*Minv( k, m ) * ( diff( M( m, i ), q_arr{ j } ) + diff( M( m, j ), q_arr{ i } ) - diff( M( i, j ), q_arr{ m } ) );
           end
           Gamma( i, j, k ) = simplify( tmp );
        end
    end
end

%% (1B) The Stiffness Terms 
syms k11 k12 k22 
syms q01 q02

Phi = 1/2* [ q01-q1, q02-q2 ]* [ k11, k12; k12, k22 ] *  [ q01-q1; q02-q2 ];
Phi = subs( Phi, { q1, q2 }, { 'q1', 'q2' } );
Phi = Phi( t );
% The stiffness components
K = sym( zeros( 2, 2 ) );

for i = 1 : 2
    for j = 1 : 2
        tmp = 0;
        for k = 1 : 2
            tmp = tmp - Gamma( i, j, k )*diff( Phi, q_arr{ k } );
        end
        K( i, j ) = simplify( diff( diff( Phi, q_arr{ i } ), q_arr{ j } ) - tmp );
    end
end


%% (1C) Initialization for the Simulation

% Substitute the Equations
old_var = { diff( q1, t, t ), diff( q2, t, t ), diff( q1, t ), diff( q2, t ),   q1,   q2 };
new_var = {           'ddq1',           'ddq2',         'dq1',         'dq2', 'q1', 'q2' };

EOM1 = subs( eq1, old_var, new_var );
EOM2 = subs( eq2, old_var, new_var );

% Getting the Mass matrix and other terms 
syms ddq1 ddq2 
b = simplify( [ EOM1; EOM2 ] - M * [ ddq1; ddq2 ] );

% Substitute the Inertia/Geometric Parameters
old_var = { M1, M2, L1, L2 };
new_var = {  1,  1,  1,  1 };

M_func = matlabFunction( subs( M, old_var, new_var ) );
b_func = matlabFunction( subs( b, old_var, new_var ) );

% Also for the Stiffness terms 
old_var2 = { k11, k12, k22 };
new_var2 = {  10,   0,  10 };
K_func = matlabFunction( subs( K, [ old_var, old_var2 ], [ new_var, new_var2 ] ) );

%% (1D) Main Simulation

T     = 5.;
dt    = 1e-4;
t_arr = 0:dt:T;
N     = length( t_arr );
t0i   = 1.0;

 q_init = [ 0, 0 ];
dq_init = [ 0, 0 ];

qi = q_init;
qf = [ 1, 1 ];
 D = 2.0;

 q_curr =  q_init;
dq_curr = dq_init;

for i = 1 : N
    
    t = t_arr( i );

    q1  =  q_curr( 1 );  q2 =  q_curr( 2 );
    dq1 = dq_curr( 1 ); dq2 = dq_curr( 2 );

    Mmat = M_func( q2 );
    bvec = b_func( dq1, dq2, q2 );

    % Torque input 
    if t <= t0i 
         q0_arr = qi;
        dq0_arr = zeros( 1, 2 );

    elseif t > t0i && t <= t0i + D
         q0_arr =  qi + ( qf - qi )* ( 10 * ( ( t-t0i )/D )^3 - 15 * ( ( t-t0i )/D )^4 +  6*( ( t-t0i )/D )^5 );
        dq0_arr = 1./D* ( qf - qi )* ( 30 * ( ( t-t0i )/D )^2 - 60 * ( ( t-t0i )/D )^3 + 30*( ( t-t0i )/D )^4 );

    else
         q0_arr = qf;
        dq0_arr = zeros( 1, 2 );

    end

    tau_in = 
    a = inv( M ) * ( - bvec  )


end
