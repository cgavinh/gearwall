import density_cube
import gaussian_tools as gt
import glob
if __name__ == "__main__":
    #c = density_cube.cube('../hocl/clhocl/clhocl_0_sp.cube')
    f = gt.FchkInterpreter(fchks=sorted(glob.glob('../hocl/clhocl/clhocl*sp.fchk')),
                           cubes=True)


    print("hello")
    print('hello')