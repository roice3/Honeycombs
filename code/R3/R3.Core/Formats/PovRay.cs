﻿namespace R3.Core
{
	using System;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using System.Text;
	using R3.Geometry;

	public class PovRay
	{
		public class Parameters
		{
			public double Scale = 1;
			public bool Halfspace = false;
			public bool ThinEdges = false;
			public double AngularThickness = 0.13;
		}

		/// <summary>
		/// Make a povray file for all the edges of an H3 model.
		/// Input edge locations are expected to live in the ball model.
		/// </summary>
		public static void WriteH3Edges( Parameters parameters, H3.Cell.Edge[] edges, string fileName, bool append = false )
		{
			WriteEdges( parameters, Geometry.Hyperbolic, edges, fileName, append );
		}

		/// <summary>
		/// Make a povray file for all the edges of a model.
		/// Works for all geometries (in conformal models, e.g. Ball and Stereographic).
		/// </summary>
		public static void WriteEdges( Parameters parameters, Geometry g, H3.Cell.Edge[] edges, string fileName, bool append )
		{
			if( append )
			{
				using( StreamWriter sw = File.AppendText( fileName ) )
				{
					foreach( H3.Cell.Edge edge in edges )
						sw.WriteLine( Edge( parameters, g, edge ) );
				}
			}
			else
			{
				using( StreamWriter sw = File.CreateText( fileName ) )
				{
					foreach( H3.Cell.Edge edge in edges )
						sw.WriteLine( Edge( parameters, g, edge ) );
				}
			}
		}

		private static string Edge( Parameters parameters, Geometry g, H3.Cell.Edge edge )
		{
			Vector3D v1 = edge.Start, v2 = edge.End;

			Vector3D[] points = null;
			Func<Vector3D, Sphere> sizeFunc = v => new Sphere() { Center = v, Radius = H3Models.SizeFuncConst( v, parameters.Scale ) };

			if( parameters.Halfspace )
			{
				v1 = H3Models.BallToUHS( v1 );
				v2 = H3Models.BallToUHS( v2 );

				points = H3Models.UHS.GeodesicPoints( v1, v2 );
				if( !parameters.ThinEdges )
					sizeFunc = v => 
					{
						// XXX, inexact
						return new Sphere() { Center = v, Radius = H3Models.UHS.SizeFunc( v, parameters.AngularThickness ) };
					};
			}
			else
			{
				if( g == Geometry.Hyperbolic )
					points = H3Models.Ball.GeodesicPoints( v1, v2 );
				else if( g == Geometry.Spherical )
				{
					points = S3.GeodesicPoints( v1, v2 );
					//points = points.Select( p => { p.Normalize(); return p; } ).ToArray();
				}
				else
				{
					//points = new Vector3D[] { v1, v2 };
					List<Vector3D> interpolated = new List<Vector3D>();
					int count = 20;
					for( int i = 0; i <= count; i++ )
						interpolated.Add( v1 + ( v2 - v1 ) * ( (double)i / count ) );
					points = interpolated.ToArray();
				}

				if( !parameters.ThinEdges )
					sizeFunc = v =>
					{
						Vector3D c;
						double r;
						H3Models.Ball.DupinCyclideSphere( v, parameters.AngularThickness/2, g, out c, out r );
						return new Sphere() { Center = c, Radius = r };
						//return new Sphere() { Center = v, Radius = H3Models.Ball.SizeFunc( v, parameters.AngularThickness ) }; // inexact
					};
			}

			//if( g == Geometry.Euclidean )
			//	return EdgeCylinder( points, sizeFunc );

			return EdgeSphereSweep( points, sizeFunc );
		}

		private static string EdgeCylinder( Vector3D[] points, Func<Vector3D, Sphere> sphereFunc )
		{
			if( points.Length != 2 )
				throw new System.ArgumentException();

			Vector3D v1 = points[0];
			Vector3D v2 = points[1];
			return string.Format( "cylinder {{ <{0:G6},{1:G6},{2:G6}>,<{3:G6},{4:G6},{5:G6}>,{6:G6} texture {{tex}} }}",
				v1.X, v1.Y, v1.Z, v2.X, v2.Y, v2.Z, sphereFunc( v1 ).Radius );
			
				// Sarah and Percy sitting in a tree...  K I S S I N G.  First comes love, then comes Roice, then goes Roice out the door and Sarah and kitty snuggling. 
				//  ^^
				// >..< 
		}

		public static string EdgeSphereSweep( Vector3D[] points, Func<Vector3D,Sphere> sphereFunc )
		{
			if( points.Length < 2 )
				throw new System.ArgumentException();

			// For the cubic spline, repeat the first and last points.
			List<Vector3D> appended = new List<Vector3D>();
			appended.Add( points.First() );
			appended.AddRange( points );
			appended.Add( points.Last() );
	
			Func<Vector3D,string> formatVecAndSize = v =>
			{
				Sphere s = sphereFunc( v );
				return string.Format( "<{0:G6},{1:G6},{2:G6}>,{3:G6}", s.Center.X, s.Center.Y, s.Center.Z, s.Radius );
			};

			// Use b_spline http://bugs.povray.org/task/81
			// options: b_spline, linear_spline, cubic_spline

			string formattedPoints = string.Join( ",", appended.Select( formatVecAndSize ).ToArray() );
			return string.Format( "sphere_sweep {{ b_spline {0}, {1} texture {{tex}} }}", points.Length + 2, formattedPoints );
		}

		/// <summary>
		/// Append facets to the end of an existing povray file.
		/// </summary>
		public static void AppendFacets( Sphere[] facets, string fileName )
		{
			using( StreamWriter sw = File.AppendText( fileName ) )
			{
				foreach( Sphere sphere in facets )
					sw.WriteLine( H3Facet( sphere ) );
			}
		}

		private static string H3Facet( Sphere sphere )
		{
			if( sphere.IsPlane )
			{
				Vector3D offsetOnNormal = Euclidean2D.ProjectOntoLine( sphere.Offset, new Vector3D(), sphere.Normal );
				return string.Format( "plane {{ {0}, {1:G6} material {{ sphereMat }} clipped_by {{ ball }} }}",
					FormatVec( sphere.Normal ), offsetOnNormal.Abs() );
			}
			else
			{
				return string.Format( "sphere {{ {0}, {1:G6} material {{ sphereMat }} clipped_by {{ ball }} }}",
					FormatVec( sphere.Center ), sphere.Radius );
			}
		}

		/// <summary>
		/// An alternative version for facets that require extra clipping.
		/// </summary>
		public static void AppendFacets( H3.Cell[] cells, string fileName )
		{
			HashSet<Sphere> completed = new HashSet<Sphere>();
			using( StreamWriter sw = File.AppendText( fileName ) )
			{
				foreach( H3.Cell cell in cells )
					sw.WriteLine( H3Facet( cell, completed ) );
			}
		}

		private static string H3Facet( H3.Cell cell, HashSet<Sphere> completed )
		{
			StringBuilder sb = new StringBuilder();

			foreach( H3.Cell.Facet facet in cell.Facets )
			{
				//if( completed.Contains( facet.Sphere ) )
				//	continue;

				// XXX - Hard coding for 535 skew
				string matString = "sphereMat";
				/*if( facet.Verts.Length == 6 )
					matString = "sphereMat";
				else if( facet.Verts.Length == 5 )
				{
					continue;
					//matString = "sphereMat2";
				}
				else
				{
					throw new System.ArgumentException();
					continue;
				}*/

				bool invert1 = !facet.Sphere.IsPointInside( cell.Center );
				if( facet.Sphere.Invert ) invert1 = !invert1;
				//bool invert1 = CheckForInvert( facet.Sphere, cell.Center );
				sb.Append( string.Format( "{0} material {{ sphereMat }} clipped_by {{ ball }}",
					FormatSphereNoMaterial( facet.Sphere, invert1, false ) ) );

				H3.Cell.Facet[] others = cell.Facets.Except( new H3.Cell.Facet[] { facet } ).ToArray();
				foreach( H3.Cell.Facet otherFacet in others )
				{
					bool invert = !otherFacet.Sphere.IsPointInside( cell.Center );
					if( otherFacet.Sphere.Invert ) invert = !invert;
					//bool invert = CheckForInvert( otherFacet.Sphere, cell.Center );
					sb.Append( string.Format( " clipped_by {{ {0} }}", FormatSphereNoMaterial( otherFacet.Sphere, invert ) ) );
				}

				sb.AppendLine( " }" );

				completed.Add( facet.Sphere );
			}

			return sb.ToString();
		}

		/// <summary>
		/// A version for the fundamental simplex.
		/// </summary>
		public static void CreateSimplex( Sphere[] facets, string fileName )
		{
			using( StreamWriter sw = File.CreateText( fileName ) )
			{
				Vector3D dummy = new Vector3D();
				sw.WriteLine( SimplexFacets( facets, dummy, new int[] { 0, 1, 2, 3 } ) );
			}
		}

		public static void AppendSimplex( Sphere[] facets, Vector3D interiorPoint, int[] include, string fileName )
		{
			using( StreamWriter sw = File.AppendText( fileName ) )
			{
				sw.WriteLine( SimplexFacets( facets, interiorPoint, include ) );
			}
		}

		private static bool CheckForInvert( Sphere facet, Vector3D point )
		{
			System.Diagnostics.Debug.Assert( !facet.IsPointInside( point ) );
			if( !facet.IsPlane )
			{
				bool invert = !facet.Invert;

				// This check doesn't seem like it should be here.
				//if( facet.IsPointInside( point ) )
				//	invert = !invert;

				return invert;
			}
			else
			{
				return !facet.Invert;
			}
		}

		private static string SimplexFacets( Sphere[] facets, Vector3D interiorPoint, int[] include )
		{
			StringBuilder sb = new StringBuilder();

			foreach( int idx in include )
			{
				Sphere facet = facets[idx];

				bool invert = CheckForInvert( facet, interiorPoint );
				sb.Append( string.Format( "{0} material {{ sphereMat2 }} clipped_by {{ ball }}",
					FormatSphereNoMaterial( facet, invert, false ) ) );

				Sphere[] others = facets.Except( new Sphere[] { facet } ).ToArray();
				foreach( Sphere otherFacet in others )
				{
					invert = CheckForInvert( otherFacet, interiorPoint );
					sb.Append( string.Format( " clipped_by {{ {0} }}", FormatSphereNoMaterial( otherFacet, invert ) ) );
				}

				sb.Append( " }" );
			}

			return sb.ToString();
		}

		private static string FormatSphereNoMaterial( Sphere sphere, bool invert, bool includeClosingBracket = true )
		{
			if( sphere.IsPlane )
			{
				Vector3D offsetOnNormal = Euclidean2D.ProjectOntoLine( sphere.Offset, new Vector3D(), sphere.Normal );
				double offset = offsetOnNormal.Abs();
				if( offsetOnNormal.Dot( sphere.Normal ) < 0 )
					offset *= -1;
				return string.Format( "plane {{ {0}, {1:G6}{2} {3}",
					FormatVec( sphere.Normal ), offset, invert ? " inverse" : string.Empty, 
					includeClosingBracket ? "}" : string.Empty );
			}
			else
			{
				return string.Format( "sphere {{ {0}, {1:G6}{2} {3}",
					FormatVec( sphere.Center ), sphere.Radius, invert ? " inverse" : string.Empty,
					includeClosingBracket ? "}" : string.Empty );
			}
		}

		public static string Sphere( Sphere sphere )
		{
			//return string.Format( "sphere {{ {0}, {1:G6} material {{ sphereMat2 }} }}",
			//		FormatVec( sphere.Center ), sphere.Radius );
			return string.Format( "sphere {{ {0}, rad material {{ sphereMat }} }}",
				FormatVec( sphere.Center ) );
		}

		private static string FormatVec( Vector3D v )
		{
			return string.Format( "<{0:G6},{1:G6},{2:G6}>", v.X, v.Y, v.Z );
		}
	}
}