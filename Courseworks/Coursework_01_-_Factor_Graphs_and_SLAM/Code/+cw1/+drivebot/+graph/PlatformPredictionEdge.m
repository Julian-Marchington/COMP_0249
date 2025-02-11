classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIAL ESTIMATE
            %
            % Provide an initial guess for x_{k+1} given x_{k}.
            % Since this factor has a measurement of zero (the mean of 
            % the process noise), we can simply set x_{k+1} = x_{k} 
            % (or include a motion model if we had one).

            % Retrieve the current state (vertex slot 1) 
            xk = obj.edgeVertices{1}.x;   % xk = [x; y; theta]
            
            % Retrieve the control input (odometry measurement)
            % u is a 3x1 vector: [v_x; v_y; omega]
            u = obj.measurement;  
            
            % Extract the current heading (theta) from xk
            theta = xk(3);
            
            % Construct the rotation matrix based on the current heading
            M = [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
            
            % Compute the predicted state at time k+1:
            % Multiply the control input by the time step dT to get the displacement.
            xk1_pred = xk + M * (u * obj.dT);
            
            % Set the initial estimate of the second vertex (state at time k+1)
            obj.edgeVertices{2}.x = xk1_pred;
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat
                       
            % Retrieve the states for time k and k+1
            xk  = obj.edgeVertices{1}.x;  % state at time k: [x; y; theta]
            xk1 = obj.edgeVertices{2}.x;  % state at time k+1
            
            % Retrieve the control input (odometry measurement)
            u = obj.measurement;   % [v_x; v_y; omega]
            
            % Get the heading from the state at time k (for constructing M)
            theta = xk(3);
            
            % Construct the rotation matrix M based on the current heading
            M = [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
            
            % Compute the error. In the absence of noise, we expect:
            %   inv(M) * (xk1 - xk) = u * dT.
            % Hence, the error is defined as:
            error = inv(M) * (xk1 - xk) - (u * obj.dT);
            
            % Store the computed error in the edge
            obj.errorZ = error;
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            obj.J{1} = -eye(3);  % derivative with respect to the first vertex (time k)
            obj.J{2} =  eye(3);  % derivative with respect to the second vertex (time k+1)
        end
    end    
end