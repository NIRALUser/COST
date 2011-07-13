#ifndef DEF_SPHERE_IKOSAHEDRON_IMAGE
#define DEF_SPHERE_IKOSAHEDRON_IMAGE

#include <itkDataObject.h>
#include <itkObjectFactory.h>
#include <itkVectorImage.h>

#include "SphereIkosahedron.h"

namespace itk {

template <class TPixel=double> 
		class ITK_EXPORT SphereIkosahedronImage : public VectorImage<TPixel,3>
{
	public:
		typedef SphereIkosahedronImage<double>			Self;
		typedef VectorImage<TPixel,3> 				Superclass;
		typedef SmartPointer<Self> 				Pointer;
		typedef SmartPointer<const Self> 			ConstPointer;
		typedef SphereIkosahedron<double>			SphereIkosahedronType;
		
		itkNewMacro(Self);
		
		void SetSubdivisionLevel(unsigned short level) 	{
			m_SphereIkosahedronObject->SetSubdivisionLevel(level);
		}
		unsigned GetSubdivisionLevel() 			{
			return m_SphereIkosahedronObject->GetSubdivisionLevel();
		}

	protected:
		SphereIkosahedronType::Pointer	m_SphereIkosahedronObject;
		SphereIkosahedronImage() 	{ m_SphereIkosahedronObject = SphereIkosahedronType::New(); };
		virtual 			~SphereIkosahedronImage() {}
		
	private:
		SphereIkosahedronImage(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
};
} // end namespace itk

#define ITK_TEMPLATE_VectorImage(_, EXPORT, x, y) namespace itk { \
 _(2(class EXPORT VectorImage< ITK_TEMPLATE_2 x >)) \
namespace Templates { typedef VectorImage< ITK_TEMPLATE_2 x > VectorImage##y; } \
}

#if ITK_TEMPLATE_TXX
# include "itkVectorImage.txx"
#endif

#ifndef ITK_MANUAL_INSTANTIATION
#include "SphereIkosahedronImage.txx"
#endif

#endif
