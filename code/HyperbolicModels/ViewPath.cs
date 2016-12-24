namespace R3.Geometry
{
	using MathNet.Numerics.Interpolation;
	using R3.Geometry;
	using System.Collections.Generic;
	using System.Linq;

	public class ViewPath
	{
		public double Time { get; set; }

		public int Step { get; set; }

		public Vector3D Location 
		{ 
			get
			{
				return new Vector3D(
					locX.Interpolate( Time ),
					locY.Interpolate( Time ),
					locZ.Interpolate( Time ) );			
			}
		}

		public Vector3D LookAt
		{ 
			get
			{
				Vector3D deriv = new Vector3D(
					locX.Differentiate( Time ), 
					locY.Differentiate( Time ), 
					locZ.Differentiate( Time ) );
				deriv.Normalize();
				return deriv;

				/*return new Vector3D(
					lookX.Interpolate( Time ),
					lookY.Interpolate( Time ),
					lookZ.Interpolate( Time ) );	*/		
			}
		}

		CubicSpline locX, locY, locZ;
		//CubicSpline lookX, lookY, lookZ;

		/// <summary>
		/// Initialize our path with a sequence of locations and lookAt directions.
		/// </summary>
		public void Initialize( IEnumerable<Vector3D> locations, IEnumerable<Vector3D> lookAts )
		{
			int count = locations.Count();
			//if( lookAts.Count() != count )
			//	throw new System.ArgumentException();
			double[] times = Enumerable.Range( 0, count ).Select( i=> (double)i/(count-1) ).ToArray();

			//.6, 0, -.8
			locX = CubicSpline.InterpolateBoundaries( times, locations.Select( v => v.X ), SplineBoundaryCondition.FirstDerivative, .6, SplineBoundaryCondition.FirstDerivative, .6 );
			locY = CubicSpline.InterpolateBoundaries( times, locations.Select( v => v.Y ), SplineBoundaryCondition.FirstDerivative, 0, SplineBoundaryCondition.FirstDerivative, 0 );
			locZ = CubicSpline.InterpolateBoundaries( times, locations.Select( v => v.Z ), SplineBoundaryCondition.FirstDerivative, -.8, SplineBoundaryCondition.FirstDerivative, -.8 );

			locX = CubicSpline.InterpolateNatural( times, locations.Select( v => v.X ) );
			locY = CubicSpline.InterpolateNatural( times, locations.Select( v => v.Y ) );
			locZ = CubicSpline.InterpolateNatural( times, locations.Select( v => v.Z ) );
			/*lookX = CubicSpline.InterpolateNatural( times, lookAts.Select( v => v.X ) );
			lookY = CubicSpline.InterpolateNatural( times, lookAts.Select( v => v.Y ) );
			lookZ = CubicSpline.InterpolateNatural( times, lookAts.Select( v => v.Z ) );*/
		}
	}
}
