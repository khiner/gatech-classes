#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix T and modifies it 
# in place to Hessenberg form using Householder reduction.
function hessenberg_form!(T)
    m = size(T, 1)

    for k = 1:m-2
        x = T[(k+1):m, k]
        v = x + sign(x[1]) * norm(x) * I[:, 1]
        v = v / norm(v)
        
        T[(k+1):m, k:m] -= 2 * v * (v' * T[(k+1):m, k:m])
        T[:, k+1:m] -= 2 * (T[:, k+1:m] * v) * v'
    end
end

function givens(a, b)
    r = sqrt(a^2 + b^2)
    c = a / r
    s = -b / r
    
    return c, s
end

#----------------------------------------
# Problem b
#----------------------------------------
# This funciton takes in a matrix T in Hessenberg form
# and runs a single iteration of the unshifted QR Algorithm 
# using Givens rotations
function givens_qr!(T)
    m = size(T, 1)
    
    for k = 1:m-1
        c, s = givens(T[k, k], T[k+1, k])
        
        G = I(m)
        G[k:k+1, k:k+1] = [c s; -s c]
        
        T = G * T
        T = T * G'
    end
    
    return T
end
    

#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix T in Hessenberg form and 
# implements the practical QR algorithm with shifts. 
# The input shift dictates which shift type your 
# algorithm should use. For shift = "single" implement the single shift 
# and for shift = "wilkinson" implement the Wilkinson shift
function practical_QR_with_shifts!(T, shift)
    max_iter = 100
    tol = 1e-6
    m = size(T, 1)

    for iter = 1:max_iter
        if shift == "single"
            mu = T[m, m]
        elseif shift == "wilkinson"
            delta = (T[m-1, m-1] - T[m, m]) / 2
            mu = T[m, m] - sign(delta) * T[m, m-1]^2 / (abs(delta) + sqrt(delta^2 + T[m, m-1]^2))
        else
            error("Invalid shift type")
        end

        T -= mu * I(m)
        T = givens_qr!(T)
        T += mu * I(m)
        
        for i = 1:m-1
            if abs(T[i+1, i]) < tol * (abs(T[i, i]) + abs(T[i+1, i+1]))
                T[i+1, i] = 0
            end
        end

        if all(T[2:end, 1:end-1] .≈ 0)
            break
        end
    end

    return T
end

#----------------------------------------
# Problem d
#----------------------------------------


