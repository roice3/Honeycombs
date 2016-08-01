namespace HyperbolicModels
{
	using R3.Core;
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.Linq;
	using System.Numerics;

	using Math = System.Math;

	public class H3Ruled
	{
		public void GenPovRay()
		{
			//H3.Cell.Edge[] fibers = Helicoid();
			H3.Cell.Edge[] fibers = Hyperboloid();
			PovRay.WriteH3Edges( new PovRay.Parameters { AngularThickness = 0.015 }, fibers, "ruled.pov" );
		}

		public H3.Cell.Edge[] Hyperboloid()
		{
			// Draw to circles of fibers, then twist them.
			List<H3.Cell.Edge> fiberList = new List<H3.Cell.Edge>();

			Vector3D cen = new Vector3D( 0, 0, 0.5 );
			double rad = .3;
			Circle3D c1 = new Circle3D { Center = cen, Radius = rad };
			Circle3D c2 = new Circle3D { Center = -cen, Radius = rad };

			int n = 50;
			Vector3D[] points1 = c1.Subdivide( n );
			Vector3D[] points2 = c2.Subdivide( n );

			double twist = 2 * Math.PI / 3;
			for( int i = 0; i < points2.Length; i++ )
			{
				points2[i].RotateXY( twist );

				Vector3D e1, e2;
				H3Models.Ball.GeodesicIdealEndpoints( points1[i], points2[i], out e1, out e2 );

				e1 = Transform( e1 );
				e2 = Transform( e2 );

				fiberList.Add( new H3.Cell.Edge( e1, e2 ) );
			}

			return fiberList.ToArray();
		}

		public H3.Cell.Edge[] Helicoid()
		{
			List<H3.Cell.Edge> fiberList = new List<H3.Cell.Edge>();

			double rotationRate = Math.PI / 40;
			int numFibers = 350;

			// Note: we need to increment a constant hyperbolic distance each step.
			int count = 0;
			double max = DonHatch.e2hNorm( 0.99 );
			double offset = max * 2 / (numFibers - 1);
			for( double z_h = -max; z_h <= max; z_h += offset )
			{
				double z = DonHatch.h2eNorm( z_h );

				Sphere s = H3Models.Ball.OrthogonalSphereInterior( new Vector3D( 0, 0, z ) );
				Circle3D c = H3Models.Ball.IdealCircle( s );

				// Two endpoints of our fiber.
				Vector3D v1 = new Vector3D( c.Radius, 0, c.Center.Z );
				Vector3D v2 = new Vector3D( -c.Radius, 0, c.Center.Z );

				v1.RotateXY( rotationRate * count );
				v2.RotateXY( rotationRate * count );

				v1 = Transform( v1 );
				v2 = Transform( v2 );

				fiberList.Add( new H3.Cell.Edge( v1, v2 ) );
				count++;
			}

			return fiberList.ToArray();
		}

		public Vector3D Transform( Vector3D v )
		{
			v.RotateAboutAxis( new Vector3D( 1, 0 ), Math.PI / 2 );

			Mobius m = new Mobius();
			m.Isometry( Geometry.Hyperbolic, 0, new Complex( 0, 0.6 ) );
			v = H3Models.TransformHelper( v, m );

			v.RotateAboutAxis( new Vector3D( 1, 0 ), -Math.PI / 2 );
			return v;
		}
	}
}
