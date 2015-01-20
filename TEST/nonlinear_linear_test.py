# -*- python -*-
# $RCSfile: nonlinear_linear_test.py,v $
# $Revision: 1.6 $
# $Author: langer $
# $Date: 2011/04/14 21:01:43 $

# This software was produced by NIST, an agency of the U.S. government,
# and by statute is not subject to copyright in the United States.
# Recipients of this software assume all responsibilities associated
# with its operation, modification and maintenance. However, to
# facilitate maintenance we ask that before distributing modified
# versions of this software, you first contact the authors at
# oof_manager@nist.gov. 


# Tests of the nonlinear solvers on linear problems.

import unittest, os
import memorycheck
import math
from UTILS import file_utils
# file_utils.generate = True

class OOF_LinearDiffusion(unittest.TestCase):
    # Solve a linear diffusion problem with a nonlinear solver.  The
    # reference output files were generated by a linear solver.

    def setUp(self):
        OOF.Microstructure.New(
            name='microstructure',
            width=1.0, height=1.0, width_in_pixels=10, height_in_pixels=10)
        OOF.Material.New(
            name='material', material_type='bulk')
        OOF.Material.Add_property(
            name='material', property='Thermal:Conductivity:Isotropic')
        OOF.Material.Add_property(
            name='material',
            property='Thermal:HeatCapacity:ConstantHeatCapacity')
        OOF.Material.Assign(
            material='material', microstructure='microstructure', pixels=all)
        OOF.Skeleton.New(
            name='skeleton', microstructure='microstructure', 
            x_elements=8, y_elements=8,
            skeleton_geometry=QuadSkeleton(
                left_right_periodicity=False,top_bottom_periodicity=False))
        OOF.Mesh.New(
            name='mesh',
            skeleton='microstructure:skeleton',
            element_types=['D2_2', 'T3_3', 'Q4_4'])
        OOF.Subproblem.Field.Define(
            subproblem='microstructure:skeleton:mesh:default',
            field=Temperature)
        OOF.Subproblem.Field.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            field=Temperature)
        OOF.Mesh.Field.In_Plane(
            mesh='microstructure:skeleton:mesh', field=Temperature)
        OOF.Subproblem.Equation.Activate(
            subproblem='microstructure:skeleton:mesh:default', 
            equation=Heat_Eqn)
        OOF.Mesh.Boundary_Conditions.New(
            name='bc',
            mesh='microstructure:skeleton:mesh',
            condition=DirichletBC(
                field=Temperature,field_component='',
                equation=Heat_Eqn,eqn_component='',
                profile=ConstantProfile(value=0.0),
                boundary='bottomleft'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<2>',
            mesh='microstructure:skeleton:mesh',
            condition=DirichletBC(
                field=Temperature,field_component='',
                equation=Heat_Eqn,eqn_component='',
                profile=ConstantProfile(value=1),
                boundary='topright'))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh', 
            field=Temperature,
            initializer=ConstScalarFieldInit(value=0.0))
        OOF.Mesh.Apply_Field_Initializers_at_Time(
            mesh='microstructure:skeleton:mesh',
            time=0.0)
        OOF.Mesh.Apply_Field_Initializers_at_Time(
            mesh='microstructure:skeleton:mesh',
            time=0.0)

        OOF.Mesh.Scheduled_Output.New(
            mesh='microstructure:skeleton:mesh',
            name='temperature', 
            output=BulkAnalysis(
                output_type='Scalar',
                data=getOutput(
                    'Field:Component', component='', field=Temperature),
                operation=DirectOutput(),
                domain=EntireMesh(),
                sampling=GridSampleSet(
                    x_points=9,y_points=9,show_x=True,show_y=True)))
        OOF.Mesh.Scheduled_Output.Schedule.Set(
            mesh='microstructure:skeleton:mesh',
            output='temperature',
            scheduletype=AbsoluteOutputSchedule(),
            schedule=Periodic(delay=0.0,interval=0.1))
        OOF.Mesh.Scheduled_Output.Destination.Set(
            mesh='microstructure:skeleton:mesh',
            output='temperature', 
            destination=OutputStream(filename='tempout.dat',mode='w'))

        OOF.Mesh.Scheduled_Output.New(
            mesh='microstructure:skeleton:mesh',
            name=AutomaticName('Average Temperature on top'),
            output=BoundaryAnalysis(operation=AverageField(field=Temperature),
                                    boundary='top'))
        OOF.Mesh.Scheduled_Output.Schedule.Set(
            mesh='microstructure:skeleton:mesh',
            output=AutomaticName('Average Temperature on top'),
            scheduletype=AbsoluteOutputSchedule(),
            schedule=Periodic(delay=0.0,interval=0.1))
        OOF.Mesh.Scheduled_Output.Destination.Set(
            mesh='microstructure:skeleton:mesh',
            output=AutomaticName('Average Temperature on top'),
            destination=OutputStream(filename='test.dat',mode='w'))

    @memorycheck.check("microstructure")
    def RKlinear(self):
        # This test uses a linear solver to generate the reference
        # data for the subsequent nonlinear solvers.
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-6,
                    minstep=1.e-05,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=RK4())),
                nonlinear_solver=NoNonlinearSolver(),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-6))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-5))
        file_utils.remove('tempout.dat')

    @memorycheck.check("microstructure")
    def CN(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-7,
                    minstep=1.e-05,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=CrankNicolson())),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-6))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-5))
        file_utils.remove('tempout.dat')

    @memorycheck.check("microstructure")
    def SS22(self):
        # Setting the AdaptiveDriver tolerance to 1.e-7 allows this
        # test to pass with the fp_file_compare tolerance set to
        # 1.e-6, but then the test runs really slowly.  Using
        # tolerances of 1.e-5 and 1.e-4, respectively, lets the test
        # run in a reasonable amount of time.
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-5,
                    minstep=1.e-6,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=SS22(theta1=0.5,theta2=0.5))),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-4))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-4))
        file_utils.remove('tempout.dat')

    @memorycheck.check("microstructure")
    def RK4(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-6,
                    minstep=1.e-05,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=RK4())),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-6))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-5))
        file_utils.remove('tempout.dat')

    # The BackwardEuler and ForwardEuler tests are done at low
    # resolution and short times, because they're slow.

    @memorycheck.check("microstructure")
    def FE(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-6,
                    minstep=1.e-06,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=ForwardEuler())),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=1.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-4,
                nlines=16))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-3,
                nlines=910))
        file_utils.remove('tempout.dat')

    @memorycheck.check("microstructure")
    def BEpicard(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-5,
                    minstep=1.e-05,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=BackwardEuler())),
                nonlinear_solver=Picard(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=1.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-3,
                nlines=16))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-3,
                nlines=910))
        file_utils.remove('tempout.dat')

    @memorycheck.check("microstructure")
    def BEnewton(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    initialstep=0.1,
                    tolerance=1.e-5,
                    minstep=1.e-05,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(singlestep=BackwardEuler())),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=1.0)

        self.assert_(file_utils.fp_file_compare(
                'test.dat',
                os.path.join('mesh_data', 'nldiff.dat'),
                tolerance=1.e-3,
                nlines=16))
        file_utils.remove('test.dat')
        self.assert_(file_utils.fp_file_compare(
                'tempout.dat',
                os.path.join('mesh_data', 'tempout.dat'),
                tolerance=1.e-3,
                nlines=910))
        file_utils.remove('tempout.dat')

    def tearDown(self):
        # outputdestination is imported explicitly so that it's still
        # found even if earlier tests are commented out.
        from ooflib.engine.IO import outputdestination
        outputdestination.forgetTextOutputStreams()
        OOF.Material.Delete(name='material')


class OOF_NLPlaneStress(unittest.TestCase):

    # Solve a trivial 1x1 linear elastic plane stress problem with a
    # nonlinear solver.  The in-plane dofs are all fixed by the
    # boundary conditions, so just the out-of-plane fields are being
    # solved.
    @memorycheck.check("microstructure")
    def Small(self):
        OOF.Microstructure.New(
            name='microstructure',
            width=1.0, height=1.0,
            width_in_pixels=10, height_in_pixels=10)
        OOF.Material.New(
            name='material', material_type='bulk')
        # Reset the default parameter values for isotropic elasticity.
        # This shouldn't be necessary if earlier tests clean up after
        # themselves properly.
        OOF.Property.Parametrize.Mechanical.Elasticity.Isotropic(
            cijkl=IsotropicRank4TensorCij(c11=1.0,c12=0.5))
        OOF.Material.Add_property(
            name='material',
            property='Mechanical:Elasticity:Isotropic')
        OOF.Material.Assign(
            material='material',
            microstructure='microstructure',
            pixels=all)
        OOF.Skeleton.New(
            name='skeleton',
            microstructure='microstructure',
            x_elements=1, y_elements=1,
            skeleton_geometry=QuadSkeleton(
                left_right_periodicity=False,top_bottom_periodicity=False))
        OOF.Mesh.New(
            name='mesh',
            skeleton='microstructure:skeleton',
            element_types=['D2_2', 'T3_3', 'Q4_4'])
        OOF.Subproblem.Field.Define(
            subproblem='microstructure:skeleton:mesh:default',
            field=Displacement)
        OOF.Subproblem.Field.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            field=Displacement)
        OOF.Subproblem.Equation.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            equation=Force_Balance)
        OOF.Subproblem.Equation.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            equation=Plane_Stress)
        OOF.Mesh.Boundary_Conditions.New(
            name='bc',
            mesh='microstructure:skeleton:mesh', 
            condition=DirichletBC(
                field=Displacement,field_component='x',
                equation=Force_Balance,eqn_component='x',
                profile=ConstantProfile(value=0.0),
                boundary='left'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<2>', 
            mesh='microstructure:skeleton:mesh', 
            condition=DirichletBC(
                field=Displacement,field_component='x',
                equation=Force_Balance,eqn_component='x',
                profile=ConstantProfile(value=0.1),
                boundary='right'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<3>', 
            mesh='microstructure:skeleton:mesh',
            condition=DirichletBC(
                field=Displacement,field_component='y',
                equation=Force_Balance,eqn_component='y',
                profile=ConstantProfile(value=0),
                boundary='left'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<4>', 
            mesh='microstructure:skeleton:mesh',
            condition=DirichletBC(
                field=Displacement,field_component='y',
                equation=Force_Balance,eqn_component='y',
                profile=ConstantProfile(value=0),
                boundary='right'))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh',
            field=Displacement, 
            initializer=ConstTwoVectorFieldInit(cx=0.0,cy=0.0))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh',
            field=Displacement_z,
            initializer=ConstThreeVectorFieldInit(cx=0.0,cy=0.0,cz=0.0))
        OOF.Mesh.Apply_Field_Initializers_at_Time(
            mesh='microstructure:skeleton:mesh',
            time=0.0)
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=StaticDriver(),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1.e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000),
                asymmetric_solver=BiConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=0.0)
        # Check that the average out-of-plane stress is 0.0.
        OOF.Mesh.Analyze.Average(
            mesh='microstructure:skeleton:mesh',
            time=latest,
            data=getOutput('Flux:Value',flux=Stress),
            domain=EntireMesh(),
            sampling=ElementSampleSet(order=automatic),
            destination=OutputStream(filename='test.dat', mode='w'))
        self.assert_(file_utils.compare_last(
                'test.dat',
                (0.0, 0.075, 0.025, 0.0, 0.0, 0.0, 0.0)))
        file_utils.remove('test.dat')

        # Check that the out-of-plane displacement derivative at the
        # nodes is correct.  u_zz should be -0.05 and the other
        # components should be 0.0.
        from ooflib.engine import mesh
        from ooflib.SWIG.engine import field
        msh = mesh.meshes["microstructure:skeleton:mesh"]
        msh_obj = msh.getObject()
        dispz = field.getField("Displacement_z")
        for node in msh_obj.funcnode_iterator():
            vals = [dispz.value(msh_obj, node, i) for i in range(3)]
            self.assert_(vals[0] == 0.0 and vals[1] == 0 and
                         math.fabs(vals[2]+0.05) < 1.e-13)

    def tearDown(self):
        # outputdestination is imported explicitly so that it's still
        # found even if earlier tests are commented out.
        from ooflib.engine.IO import outputdestination
        outputdestination.forgetTextOutputStreams()
        OOF.Material.Delete(name='material')


# Plane stress elasticity with a floating boundary condition.  A
# thermal diffusion version of this problem is done in
# nonlinear_floatbc.py.  Perhaps this test should be moved to that
# file.

class OOF_NLPlaneStress2(unittest.TestCase):
    def setUp(self):
        OOF.Microstructure.New(
            name='microstructure',
            width=1.0, height=1.0,
            width_in_pixels=10, height_in_pixels=10)
        OOF.Material.New(
            name='material', material_type='bulk')
        OOF.Material.Assign(
            material='material', microstructure='microstructure', pixels=all)
        # Reset the default parameter values for isotropic elasticity.
        # This shouldn't be necessary if earlier tests clean up after
        # themselves properly.
        OOF.Property.Parametrize.Mechanical.Elasticity.Isotropic(
            cijkl=IsotropicRank4TensorCij(c11=1.0,c12=0.5))
        OOF.Material.Add_property(
            name='material', property='Mechanical:Elasticity:Isotropic')
        OOF.Material.Add_property(
            name='material',
            property='Mechanical:MassDensity:ConstantMassDensity')
        OOF.Skeleton.New(
            name='skeleton',
            microstructure='microstructure',
            x_elements=5, y_elements=5,
            skeleton_geometry=QuadSkeleton(
                left_right_periodicity=False,top_bottom_periodicity=False))
        OOF.Mesh.New(
            name='mesh',
            skeleton='microstructure:skeleton',
            element_types=['D2_2', 'T3_3', 'Q4_4'])
        OOF.Subproblem.Field.Define(
            subproblem='microstructure:skeleton:mesh:default',
            field=Displacement)
        OOF.Subproblem.Field.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            field=Displacement)
        OOF.Subproblem.Equation.Activate(
            subproblem='microstructure:skeleton:mesh:default',
            equation=Force_Balance)
        OOF.Subproblem.Equation.Activate(
            subproblem='microstructure:skeleton:mesh:default', 
            equation=Plane_Stress)
        OOF.Mesh.Boundary_Conditions.New(
            name='bc',
            mesh='microstructure:skeleton:mesh',
            condition=DirichletBC(
                field=Displacement, field_component='x',
                equation=Force_Balance, eqn_component='x',
                profile=ConstantProfile(value=0.0),
                boundary='left'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<2>',
            mesh='microstructure:skeleton:mesh', 
            condition=DirichletBC(
                field=Displacement,field_component='y',
                equation=Force_Balance,eqn_component='y',
                profile=ConstantProfile(value=0.0),
                boundary='bottomleft'))
        OOF.Mesh.Boundary_Conditions.New(
            name='bc<3>', 
            mesh='microstructure:skeleton:mesh',
            condition=FloatBC(
                field=Displacement,field_component='x',
                equation=Force_Balance,eqn_component='x',
                profile=ConstantProfile(value=0.0),
                boundary='right'))
        OOF.Mesh.Scheduled_Output.New(
            mesh='microstructure:skeleton:mesh',
            name='riiight',
            output=BoundaryAnalysis(
                operation=AverageField(field=Displacement),
                boundary='right'))
        OOF.Mesh.Scheduled_Output.Schedule.Set(
            mesh='microstructure:skeleton:mesh',
            output='riiight',
            scheduletype=AbsoluteOutputSchedule(),
            schedule=Periodic(delay=0.0,interval=0.1))
        OOF.Mesh.Scheduled_Output.Destination.Set(
            mesh='microstructure:skeleton:mesh',
            output='riiight',
            destination=OutputStream(filename='riiight.out', mode='w'))
        OOF.Mesh.Scheduled_Output.New(
            mesh='microstructure:skeleton:mesh',
            name='avgstress', 
            output=BulkAnalysis(
                output_type='Aggregate',
                data=getOutput('Flux:Value',flux=Stress),
                operation=AverageOutput(),
                domain=EntireMesh(),
                sampling=ElementSampleSet(order=automatic)))
        OOF.Mesh.Scheduled_Output.Schedule.Set(
            mesh='microstructure:skeleton:mesh',
            output='avgstress', 
            scheduletype=AbsoluteOutputSchedule(),
            schedule=Periodic(delay=0.0,interval=0.1))
        OOF.Mesh.Scheduled_Output.Destination.Set(
            mesh='microstructure:skeleton:mesh',
            output='avgstress', 
            destination=OutputStream(filename='stress_nl2.out',mode='w'))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh',
            field=Displacement,
            initializer=FuncTwoVectorFieldInit(fx='0.1*x',fy='0.0'))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh',
            field=Displacement_t, 
            initializer=ConstTwoVectorFieldInit(cx=0.0,cy=0.0))
        OOF.Mesh.Set_Field_Initializer(
            mesh='microstructure:skeleton:mesh',
            field=Displacement_z, 
            initializer=ConstThreeVectorFieldInit(cx=0.0,cy=0.0,cz=0.0))
        OOF.Mesh.Apply_Field_Initializers_at_Time(
            mesh='microstructure:skeleton:mesh',
            time=0.0)
    def check(self):
        self.assert_(file_utils.fp_file_compare(
                'riiight.out',
                os.path.join('mesh_data', 'riiight.out'),
                tolerance=1.e-4))
        file_utils.remove('riiight.out')
        self.assert_(file_utils.fp_file_compare(
                'stress_nl2.out',
                os.path.join('mesh_data', 'stress_nl2.out'),
                tolerance=1.e-4))
        file_utils.remove('stress_nl2.out')

    @memorycheck.check("microstructure")
    def LinearCN(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    tolerance=1e-6,
                    initialstep=0,
                    minstep=1e-08,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(
                        singlestep=CrankNicolson())),
                nonlinear_solver=NoNonlinearSolver(),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000),
                asymmetric_solver=BiConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)
        self.check()
    @memorycheck.check("microstructure")
    def LinearSS22(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    tolerance=1e-6,
                    initialstep=0,
                    minstep=1e-08,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(
                        singlestep=SS22(theta1=0.5,theta2=0.5))),
                nonlinear_solver=NoNonlinearSolver(),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000),
                asymmetric_solver=BiConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)
        self.check()
    @memorycheck.check("microstructure")
    def NonlinearSS22(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    tolerance=1e-6,
                    initialstep=0,
                    minstep=1e-08,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(
                        singlestep=SS22(theta1=0.5,theta2=0.5))),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000),
                asymmetric_solver=BiConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)
        self.check()
    @memorycheck.check("microstructure")
    def NonlinearCN(self):
        OOF.Subproblem.Set_Solver(
            subproblem='microstructure:skeleton:mesh:default',
            solver_mode=AdvancedSolverMode(
                time_stepper=AdaptiveDriver(
                    tolerance=1e-6,
                    initialstep=0,
                    minstep=1e-08,
                    errorscaling=AbsoluteErrorScaling(),
                    stepper=TwoStep(
                        singlestep=CrankNicolson())),
                nonlinear_solver=Newton(
                    relative_tolerance=1e-08,
                    absolute_tolerance=1e-13,
                    maximum_iterations=200),
                symmetric_solver=ConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000),
                asymmetric_solver=BiConjugateGradient(
                    preconditioner=ILUPreconditioner(),
                    tolerance=1e-13,
                    max_iterations=1000)))
        OOF.Mesh.Solve(
            mesh='microstructure:skeleton:mesh',
            endtime=2.0)
        self.check()

    def tearDown(self):
        from ooflib.engine.IO import outputdestination
        outputdestination.forgetTextOutputStreams()
        OOF.Material.Delete(name='material')


#=--=##=--=##=--=##=--=##=--=##=--=##=--=##=--=##=--=##=--=##=--=#

# Routine to do regression-type testing on the items in this file.
# Tests will be run in the order they appear in the list.  This
# routine will stop after the first failure.

def run_tests():
    test_set = [
        ## RKLinear just generates the reference files for the rest of
        ## the LinearDiffusion tests.  It should come first.
        OOF_LinearDiffusion("RKlinear"),
        OOF_LinearDiffusion("CN"),
        OOF_LinearDiffusion("SS22"),
        OOF_LinearDiffusion("RK4"),
        OOF_LinearDiffusion("FE"),
        OOF_LinearDiffusion("BEpicard"),
        OOF_LinearDiffusion("BEnewton"),

        OOF_NLPlaneStress("Small"),
        OOF_NLPlaneStress2("LinearCN"),
        OOF_NLPlaneStress2("LinearSS22"),
        OOF_NLPlaneStress2("NonlinearCN"),
        OOF_NLPlaneStress2("NonlinearSS22")
        ]
    
    #test_set = [OOF_NLPlaneStress2("LinearCN")]
    # test_set = [OOF_NLPlaneStress("Small")]

    logan = unittest.TextTestRunner()
    for t in test_set:
        print >> sys.stderr,  "\n *** Running test: %s\n" % t.id()
        res = logan.run(t)
        if not res.wasSuccessful():
            return 0
    return 1


###################################################################
# The code below this line should be common to all testing files. #
###################################################################

if __name__=="__main__":
    # If directly run, then start oof, and run the local tests, then quit.
    import sys
    try:
        import oof2
        sys.path.append(os.path.dirname(oof2.__file__))
        from ooflib.common import oof
        from math import *
    except ImportError:
        print "OOF is not correctly installed on this system."
        sys.exit(4)
    sys.argv.append("--text")
    sys.argv.append("--quiet")
    sys.argv.append("--seed=17")
    oof.run(no_interp=1)

    success = run_tests()

    OOF.File.Quit()
    
    if success:
        print "All tests passed."
    else:
        print "Test failure."
