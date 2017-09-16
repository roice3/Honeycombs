namespace R3.Math
{
	using Math = System.Math;

	using R3.Geometry;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Numerics;
	using System.Runtime.Serialization;

	/// <summary>
	/// Class to represent an isometry.
	/// This is really just a wrapper around a Mobius transformation, but also includes a reflection in a generalized circle.
	/// (Reflections can't be defined with a Mobius transformation.)
	/// NOTE: The order in which the two elements are applied is important.  We will apply the Mobius part of the isometry first.
	/// </summary>
	[DataContract( Namespace = "" )]
	public class Isometry : ITransform
	{
		public Isometry()
		{
			m_mobius.Unity();
		}

		public Isometry( Mobius m, Circle r )
		{
			Mobius = m;
			Reflection = r;
		}

		public Isometry( Isometry i )
		{
			Mobius = i.Mobius;
			if( i.Reflection != null )
				Reflection = i.Reflection.Clone();
		}

		public Isometry Clone()
		{
			return new Isometry( this );
		}

		/// <summary>
		/// Mobius Transform for this isometry.
		/// </summary>
		[DataMember]
		public Mobius Mobius 
		{ 
			get { return m_mobius; } 
			set { m_mobius = value; } 
		}
		private Mobius m_mobius;

		/// <summary>
		/// Defines the circle (or line) in which to reflect for this isometry.
		/// Null if we don't want to include a reflection.
		/// </summary>
		[DataMember]
		public Circle Reflection
		{
			get { return m_reflection; }
			set 
			{ 
				m_reflection = value;
				CacheCircleInversion( m_reflection );
			}
		}

		/// <summary>
		/// Whether or not we are reflected.
		/// </summary>
		public bool Reflected
		{
			get
			{
				return m_reflection != null;
			}
		}

		// NOTE: Applying isometries with reflections was really slow, so we cache the Mobius transforms we need to more quickly do it.
		private Circle m_reflection;
		private Mobius m_cache1;
		private Mobius m_cache2;

		/// <summary>
		/// Applies an isometry to a vector.
		/// </summary>
		/// <remarks>Use the complex number version if you can.</remarks>
		public Vector3D Apply( Vector3D z )
		{
			Complex cInput = z;
			Complex cOutput = Apply( cInput );
			return Vector3D.FromComplex( cOutput );
		}

		/// <summary>
		/// Applies an isometry to a complex number.
		/// </summary>
		public Complex Apply( Complex z )
		{
			z = Mobius.Apply( z );
			if( Reflection != null )
				z = ApplyCachedCircleInversion( z );
			return z;
		}

		/// <summary>
		/// Does a circle inversion on an arbitrary circle.
		/// </summary>
		private void CacheCircleInversion( Circle inversionCircle )
		{
			if( inversionCircle == null )
				return;

			Complex p1, p2, p3;
			if( inversionCircle.IsLine )
			{
				p1 = inversionCircle.P1;
				p2 = inversionCircle.P2;
				p3 = (p1 + p2) / 2;
			}
			else
			{
				p1 = (inversionCircle.Center + new Vector3D( inversionCircle.Radius, 0 ));
				p2 = (inversionCircle.Center + new Vector3D( -inversionCircle.Radius, 0 ));
				p3 = (inversionCircle.Center + new Vector3D( 0, inversionCircle.Radius ));
			}

			CacheCircleInversion( p1, p2, p3 );
		}

		/// <summary>
		/// Does a circle inversion in an arbitrary, generalized circle.
		/// IOW, the three points may be collinear, in which case we are talking about a reflection.
		/// </summary>
		private void CacheCircleInversion( Complex c1, Complex c2, Complex c3 )
		{
			Mobius toUnitCircle = new Mobius();
			toUnitCircle.MapPoints(
				c1, c2, c3,
				new Complex( 1, 0 ),
				new Complex( -1, 0 ),
				new Complex( 0, 1 ) );

			m_cache1 = toUnitCircle;
			m_cache2 = m_cache1.Inverse();
		}

		private Complex ApplyCachedCircleInversion( Complex input )
		{
			Complex result = m_cache1.Apply( input );
			result = CircleInversion( result );
			result = m_cache2.Apply( result );
			return result;
		}

		private static bool IsNaN( Complex c )
		{
			return
				double.IsNaN( c.Real ) ||
				double.IsNaN( c.Imaginary );
		}

		/// <summary>
		/// This will reflect a point in an origin centered circle.
		/// </summary>
		private Complex CircleInversion( Complex input )
		{
			if( IsNaN( input ) )
				return Complex.Zero;

			return Complex.One / Complex.Conjugate( input );
		}
	}
}
